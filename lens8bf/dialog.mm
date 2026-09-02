// The Filmic filter dialog.
//
// Laid out the way Photoshop's own filter dialogs are, because that is what the
// muscle memory expects: preview on the left with zoom controls beneath it,
// OK/Cancel top right, Preview checkbox under them, a numeric field beside every
// slider, and Option held down turning Cancel into Reset.
//
// The preview is fetched from Photoshop and re-fetched only when the zoom or pan
// changes. Slider moves re-render from the copy already in hand: re-entering the
// host for pixels on every drag would be slow and is a re-entrancy hazard inside
// a modal loop. Rendering runs off the main thread so dragging never stutters.
//
// NOTE ON CANVAS PREVIEW: a filter plug-in cannot preview onto the document.
// Nothing in the filter API exposes it -- Photoshop's own filters can only do it
// because they are inside Photoshop. The zoomable proxy here is the equivalent.

#import <Cocoa/Cocoa.h>

#include "dialog.h"
#include "presets.h"
#include "render.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {
// Starting size of the preview well; it grows with the window from here.
constexpr int kPreviewW = 640;
constexpr int kPreviewH = 420;
constexpr CGFloat kPanelW = 330;   // right-hand control column, fixed width
constexpr CGFloat kGap    = 26;
constexpr CGFloat kMargin = 20;
// Everything below the preview: zoom row, status, preset, checkboxes, buttons.
constexpr CGFloat kBelowPreview = 210;
// Ask Photoshop for exactly the well's shape. Asking for a single "max side"
// clamps both dimensions to the same number, which crops a square out of a wide
// frame as soon as you zoom past fit.
}


// ---------------------------------------------------------------- saved presets

// Saved presets live in the plug-in's own defaults domain rather than in a file
// we invent a location for: Photoshop hosts us, so anything we write beside the
// bundle may sit inside /Applications and need admin rights.
static NSString* const kPresetsKey = @"savedPresets";

static NSUserDefaults* presetStore(void) {
    return [[NSUserDefaults alloc] initWithSuiteName:@"com.moggach.filmic"];
}

static NSDictionary* controlsToDict(const LensControls& c) {
    return @{@"distortion": @(c.distortion), @"ca": @(c.chromaticAberration),
             @"vignette": @(c.vignette),     @"edgeBlur": @(c.edgeBlur),
             @"astigmatism": @(c.astigmatism), @"coma": @(c.coma),
             @"tStop": @(c.tStop),           @"squeeze": @(c.squeeze),
             @"grain": @(c.grain),           @"grainSize": @(c.grainSize),
             @"grainColour": @(c.grainColour), @"bands": @(c.bands)};
}

static LensControls dictToControls(NSDictionary* d) {
    LensControls c{};
    auto get = [&](NSString* k, double fallback) {
        id v = d[k];
        return [v isKindOfClass:[NSNumber class]] ? [v doubleValue] : fallback;
    };
    c.distortion          = get(@"distortion", 0);
    c.chromaticAberration = get(@"ca", 0);
    c.vignette            = get(@"vignette", 0);
    c.edgeBlur            = get(@"edgeBlur", 0);
    c.astigmatism         = get(@"astigmatism", 0);
    c.coma                = get(@"coma", 0);
    c.tStop               = get(@"tStop", 2.0);
    c.squeeze             = get(@"squeeze", 1.0);
    c.grain               = get(@"grain", 0);
    c.grainSize           = get(@"grainSize", 1.5);
    c.grainColour         = get(@"grainColour", 35);
    c.bands               = int(get(@"bands", 7));
    return c;
}

// ---------------------------------------------------------------------------

// The preview well. It exists as its own class only so it can report drags: at
// anything past "fit" the view shows a window onto a larger image, and dragging
// is how you get to the corners, which is exactly where lens aberrations live.
@interface LensPreviewView : NSImageView
@property(nonatomic, copy) void (^onDrag)(CGFloat dx, CGFloat dy);
@property(nonatomic, assign) BOOL pannable;
@end

@implementation LensPreviewView {
    NSPoint _last;
}
- (void)mouseDown:(NSEvent*)e { _last = [self convertPoint:e.locationInWindow fromView:nil]; }
- (void)mouseDragged:(NSEvent*)e {
    if (!self.pannable || !self.onDrag) return;
    const NSPoint p = [self convertPoint:e.locationInWindow fromView:nil];
    self.onDrag(p.x - _last.x, p.y - _last.y);
    _last = p;
}
- (void)resetCursorRects {
    [self addCursorRect:self.bounds
                 cursor:self.pannable ? [NSCursor openHandCursor] : [NSCursor arrowCursor]];
}
@end

@interface LensDialogController : NSObject <NSWindowDelegate>
@property(nonatomic, strong) NSWindow* window;
@property(nonatomic, assign) BOOL accepted;
@end

@implementation LensDialogController {
    LensControls* _controls;
    ProxySource*  _source;

    std::vector<float> _original;      // the proxy exactly as Photoshop gave it
    int _pw, _ph;
    double _pixelScale;
    std::string _tablePath;

    int    _zoomScale;                 // image pixels per preview pixel; 1 == 100%
    int    _fitScale;
    double _centreX, _centreY;         // where we are looking, 0..1

    LensPreviewView* _preview;
    NSView*      _rightPanel;      // the whole control column, moved as one
    NSButton*    _zoomOutButton;
    NSButton*    _zoomInButton;
    NSButton*    _okButton;
    NSTextField* _presetLabel;
    int          _previewW, _previewH;
    NSSize       _lastLaidOut;
    NSTextField* _zoomLabel;
    NSButton*    _previewToggle;
    NSButton*    _gridToggle;
    NSButton*    _cancelButton;
    NSTextField* _status;
    NSMutableArray<NSSlider*>*    _sliders;
    NSMutableArray<NSTextField*>* _fields;
    NSButton* _qDraft;
    NSButton* _qNormal;
    NSButton* _qBest;
    NSBox* _lensBox;
    NSBox* _filmBox;
    NSView* _canvas;
    NSPopUpButton* _presets;
    NSTextField*   _presetDesc;
    id _modifierMonitor;

    BOOL _renderPending, _renderRunning;
    BOOL _dragging;                 // a slider is being held, so render coarse
}

- (instancetype)initWithControls:(LensControls*)c source:(ProxySource*)src {
    if ((self = [super init])) {
        _controls = c;
        _source = src;
        _centreX = _centreY = 0.5;
        _zoomScale = 0;               // 0 means "fit", resolved on first fetch
        // The first fetch happens BEFORE build lays anything out, so the well's
        // size has to be known here. Left at zero it asked Photoshop to fit the
        // picture into 64 pixels, which came back as a blocky 2% thumbnail and
        // stayed that way until the window was resized.
        _previewW = kPreviewW;
        _previewH = kPreviewH;
        _tablePath = src->tablePath();
        _accepted = NO;
        if (![self refetch]) return nil;
        [self build];
    }
    return self;
}

- (void)dealloc {
    if (_modifierMonitor) [NSEvent removeMonitor:_modifierMonitor];
}

// ------------------------------------------------------------ image fetching --

- (BOOL)refetch {
    std::vector<float> px;
    int w = 0, h = 0;
    double scale = 1.0;
    if (!_source->fetchProxy(std::max(64, _previewW), std::max(64, _previewH),
                             _zoomScale, _centreX, _centreY, px, w, h, scale))
        return NO;
    _original = std::move(px);
    _pw = w; _ph = h;
    _pixelScale = scale;
    if (_zoomScale == 0) { _zoomScale = int(scale); _fitScale = int(scale); }
    return YES;
}

// ------------------------------------------------------------------ building --

- (NSSlider*)addRowAt:(CGFloat)y x:(CGFloat)x w:(CGFloat)width
                label:(NSString*)name desc:(NSString*)desc
                  min:(double)lo max:(double)hi value:(double)v
                  tag:(NSInteger)tag into:(NSView*)host {
    // Label and value share a line, slider directly beneath, and the explanation
    // lives in a tooltip. Photoshop's own dialogs carry no standing prose under
    // each control, and eleven description lines were most of the reason this
    // column read as a different piece of software -- they roughly doubled its
    // height for text you read once.
    NSTextField* title = [NSTextField labelWithString:name];
    title.font = [NSFont systemFontOfSize:13];   // Photoshop's own UI size
    title.frame = NSMakeRect(x, y, width - 70, 17);
    title.toolTip = desc;
    [host addSubview:title];

    NSTextField* field = [[NSTextField alloc] initWithFrame:NSMakeRect(x + width - 62, y - 3, 62, 21)];
    field.alignment = NSTextAlignmentRight;
    field.font = [NSFont monospacedDigitSystemFontOfSize:12 weight:NSFontWeightRegular];
    field.target = self;
    field.action = @selector(fieldEdited:);
    field.tag = tag;
    field.toolTip = desc;
    [host addSubview:field];
    [_fields addObject:field];

    NSSlider* sl = [NSSlider sliderWithValue:v minValue:lo maxValue:hi
                                      target:self action:@selector(sliderMoved:)];
    sl.frame = NSMakeRect(x, y - 27, width, 18);
    sl.tag = tag;
    sl.continuous = YES;
    sl.toolTip = desc;
    [host addSubview:sl];
    [_sliders addObject:sl];
    return sl;
}

- (void)build {
    _sliders = [NSMutableArray array];
    _fields  = [NSMutableArray array];

    const CGFloat rowStep = 49;   // measured from Photoshop's own slider rows
    const CGFloat lensBoxH = 8 * rowStep + 42;   // room under the title
    const CGFloat filmBoxH = 3 * rowStep + 42;
    const CGFloat columnH  = 42 + 30 + 34 + lensBoxH + 14 + filmBoxH;
    const CGFloat winW = kMargin + kPreviewW + kGap + kPanelW + kMargin;
    const CGFloat winH = std::max<CGFloat>(columnH + 2 * kMargin,
                                           kPreviewH + kBelowPreview + kMargin);

    self.window = [[NSWindow alloc]
        initWithContentRect:NSMakeRect(0, 0, winW, winH)
                  styleMask:NSWindowStyleMaskTitled | NSWindowStyleMaskResizable
                    backing:NSBackingStoreBuffered
                      defer:NO];
    self.window.title = @"Filmic";
    self.window.delegate = self;
    self.window.contentMinSize = NSMakeSize(kMargin * 2 + 420 + kGap + kPanelW,
                                            columnH + 2 * kMargin);
    [self.window center];
    NSView* content = self.window.contentView;

    // ---- canvas: a recessed area the picture sits inside, as Lens Blur does --
    _canvas = [[NSView alloc] initWithFrame:NSZeroRect];
    _canvas.wantsLayer = YES;
    _canvas.layer.backgroundColor = [NSColor underPageBackgroundColor].CGColor;
    _canvas.layer.borderWidth = 1;
    _canvas.layer.borderColor = [NSColor separatorColor].CGColor;
    [content addSubview:_canvas];

    _preview = [[LensPreviewView alloc] initWithFrame:NSZeroRect];
    _preview.imageScaling = NSImageScaleProportionallyUpOrDown;   // letterbox, never stretch
    __weak LensDialogController* panSelf = self;
    _preview.onDrag = ^(CGFloat dx, CGFloat dy) { [panSelf panByX:dx y:dy]; };
    [_canvas addSubview:_preview];

    // ---- zoom, bottom left under the canvas ----
    _zoomOutButton = [NSButton buttonWithTitle:@"−" target:self action:@selector(zoomOut:)];
    _zoomOutButton.bezelStyle = NSBezelStyleCircular;
    [content addSubview:_zoomOutButton];

    _zoomLabel = [NSTextField labelWithString:@"100%"];
    _zoomLabel.alignment = NSTextAlignmentCenter;
    _zoomLabel.font = [NSFont monospacedDigitSystemFontOfSize:11 weight:NSFontWeightRegular];
    [content addSubview:_zoomLabel];

    _zoomInButton = [NSButton buttonWithTitle:@"+" target:self action:@selector(zoomIn:)];
    _zoomInButton.bezelStyle = NSBezelStyleCircular;
    [content addSubview:_zoomInButton];

    _status = [NSTextField labelWithString:@""];
    _status.font = [NSFont systemFontOfSize:10];
    _status.textColor = [NSColor secondaryLabelColor];
    [content addSubview:_status];

    _presetLabel = [NSTextField labelWithString:@"Preset"];
    _presetLabel.font = [NSFont systemFontOfSize:12];
    [content addSubview:_presetLabel];

    _presets = [[NSPopUpButton alloc] initWithFrame:NSZeroRect];
    _presets.target = self;
    _presets.action = @selector(presetChosen:);
    [content addSubview:_presets];

    _presetDesc = [NSTextField labelWithString:@""];
    _presetDesc.font = [NSFont systemFontOfSize:11];
    _presetDesc.textColor = [NSColor secondaryLabelColor];
    [content addSubview:_presetDesc];

    // ---- right column: OK and Cancel at the very top ----
    _okButton = [NSButton buttonWithTitle:@"OK" target:self action:@selector(ok:)];
    _okButton.keyEquivalent = @"\r";
    [content addSubview:_okButton];

    // One action, decided at click time. Swapping the action from the modifier
    // monitor raced the click: flags can change between press and release, so the
    // button ran whichever action arrived last and Reset did nothing.
    _cancelButton = [NSButton buttonWithTitle:@"Cancel" target:self
                                       action:@selector(cancelOrReset:)];
    _cancelButton.keyEquivalent = @"\033";
    _cancelButton.toolTip = @"Hold Option to turn this into Reset.";
    [content addSubview:_cancelButton];

    _previewToggle = [NSButton checkboxWithTitle:@"Preview" target:self action:@selector(rerender:)];
    _previewToggle.state = NSControlStateValueOn;
    [content addSubview:_previewToggle];

    _gridToggle = [NSButton checkboxWithTitle:@"Guides" target:self action:@selector(rerender:)];
    _gridToggle.state = NSControlStateValueOff;
    _gridToggle.toolTip = @"Lay a grid of straight lines and dots over the picture. "
                           "Distortion bends the lines, chromatic aberration colours them, "
                           "blur spreads the dots.";
    [content addSubview:_gridToggle];

    // Quality as one row of choices rather than a menu: three states you can see
    // at once, the way Lens Blur puts Faster against More Accurate.
    _qDraft  = [NSButton radioButtonWithTitle:@"Draft"  target:self action:@selector(rerender:)];
    _qNormal = [NSButton radioButtonWithTitle:@"Normal" target:self action:@selector(rerender:)];
    _qBest   = [NSButton radioButtonWithTitle:@"Best"   target:self action:@selector(rerender:)];
    for (NSButton* b in @[_qDraft, _qNormal, _qBest]) {
        b.toolTip = @"How many wavelengths of light are simulated: 7, 11 or 15. "
                     "More gives smoother, truer colour in the fringing, and costs "
                     "time in proportion.";
        [content addSubview:b];
    }
    _qNormal.state = NSControlStateValueOn;

    // ---- grouped boxes, titled, as Lens Blur groups Iris and Noise ----
    _lensBox = [[NSBox alloc] initWithFrame:NSMakeRect(0, 0, kPanelW, lensBoxH)];
    _lensBox.title = @"Lens";
    _lensBox.titlePosition = NSAtTop;
    _lensBox.titleFont = [NSFont boldSystemFontOfSize:13];
    [content addSubview:_lensBox];

    _filmBox = [[NSBox alloc] initWithFrame:NSMakeRect(0, 0, kPanelW, filmBoxH)];
    _filmBox.title = @"Film";
    _filmBox.titlePosition = NSAtTop;
    _filmBox.titleFont = [NSFont boldSystemFontOfSize:13];
    [content addSubview:_filmBox];

    const CGFloat innerW = kPanelW - 26;
    CGFloat y = _lensBox.contentView.frame.size.height - 34;
    // Ordered the way a lens is reached for: the shape of the frame first, then
    // the character of the blur.
    [self addRowAt:y x:0 w:innerW label:@"Anamorphic squeeze"
              desc:@"1.0 is a normal lens. 2.0 gives oval bokeh and sideways streaks."
               min:1.0 max:2.0 value:_controls->squeeze tag:0 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Distortion"
              desc:@"Minus bulges outward (barrel), plus pinches in (pincushion)."
               min:-100 max:100 value:_controls->distortion tag:1 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Vignette"
              desc:@"Darkens the corners. Minus brightens them instead."
               min:-100 max:100 value:_controls->vignette tag:2 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Chromatic aberration"
              desc:@"Colour fringing on edges away from the centre."
               min:-100 max:100 value:_controls->chromaticAberration tag:3 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Aperture"
              desc:@"Sets how big the softness gets. Lower is softer."
               min:1.4 max:16 value:_controls->tStop tag:4 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Edge softness"
              desc:@"Sharp in the middle, soft at the edges (field curvature)."
               min:-100 max:100 value:_controls->edgeBlur tag:5 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Astigmatism"
              desc:@"Stretches edge detail into streaks."
               min:-100 max:100 value:_controls->astigmatism tag:6 into:_lensBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Coma"
              desc:@"Pulls corner highlights into comet tails."
               min:-100 max:100 value:_controls->coma tag:7 into:_lensBox.contentView];

    y = _filmBox.contentView.frame.size.height - 34;
    [self addRowAt:y x:0 w:innerW label:@"Grain"
              desc:@"Strongest in the mid tones, as on real film."
               min:0 max:100 value:_controls->grain tag:8 into:_filmBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Grain size"
              desc:@"Clump size in pixels, independent of image size."
               min:0.5 max:4.0 value:_controls->grainSize tag:9 into:_filmBox.contentView]; y -= rowStep;
    [self addRowAt:y x:0 w:innerW label:@"Grain colour"
              desc:@"0 is monochrome grain, 100 speckles each channel apart."
               min:0 max:100 value:_controls->grainColour tag:10 into:_filmBox.contentView];

    // The monitor only retitles the button; the click itself decides what to do.
    __weak LensDialogController* weak = self;
    _modifierMonitor = [NSEvent addLocalMonitorForEventsMatchingMask:NSEventMaskFlagsChanged
                                                             handler:^NSEvent*(NSEvent* e) {
        LensDialogController* me = weak;
        if (me) {
            const BOOL alt = (e.modifierFlags & NSEventModifierFlagOption) != 0;
            me->_cancelButton.title = alt ? @"Reset" : @"Cancel";
        }
        return e;
    }];

    // layoutViews decides the real size of the well, which is larger than the
    // constant the first fetch used, so take the picture again at that size.
    [self layoutViews];
    [self refetch];

    [self rebuildPresetMenu];
    for (const auto& pr : lensPresets())
        if (std::string(pr.name) == "Default") _presetDesc.stringValue = @(pr.desc);
    [self syncFields];
    [self updateZoomLabel];
    [self updatePannable];
    [self rerender:nil];
}

// Frames are set in one place rather than at construction, so a resize is the
// same code path as the first layout and the two cannot drift apart.
- (void)layoutViews {
    const NSSize sz = self.window.contentView.bounds.size;
    const CGFloat panelX = sz.width - kMargin - kPanelW;
    CGFloat y = sz.height - kMargin;

    // Buttons first, at the top right, above everything else.
    _okButton.frame     = NSMakeRect(panelX + kPanelW - 106, y - 32, 106, 32);
    _cancelButton.frame = NSMakeRect(panelX + kPanelW - 218, y - 32, 106, 32);
    y -= 42;

    _previewToggle.frame = NSMakeRect(panelX, y - 20, 90, 20);
    _gridToggle.frame    = NSMakeRect(panelX + 92, y - 20, 90, 20);
    y -= 30;

    _qDraft.frame  = NSMakeRect(panelX,       y - 20, 74, 20);
    _qNormal.frame = NSMakeRect(panelX + 78,  y - 20, 82, 20);
    _qBest.frame   = NSMakeRect(panelX + 164, y - 20, 74, 20);
    y -= 34;

    _lensBox.frame = NSMakeRect(panelX, y - _lensBox.frame.size.height,
                                kPanelW, _lensBox.frame.size.height);
    y -= _lensBox.frame.size.height + 14;
    _filmBox.frame = NSMakeRect(panelX, y - _filmBox.frame.size.height,
                                kPanelW, _filmBox.frame.size.height);

    // Canvas fills everything left of the column, down to the zoom row.
    const CGFloat cw = panelX - kGap - kMargin;
    const CGFloat ch = sz.height - kMargin - kBelowPreview;
    _canvas.frame = NSMakeRect(kMargin, sz.height - kMargin - ch, cw, ch);
    _preview.frame = NSMakeRect(6, 6, cw - 12, ch - 12);
    _previewW = int(cw - 12);
    _previewH = int(ch - 12);

    // Zoom at the bottom LEFT, under the canvas, where Lens Blur keeps it.
    const CGFloat zoomY = _canvas.frame.origin.y - 34;
    _zoomOutButton.frame = NSMakeRect(kMargin, zoomY, 26, 26);
    _zoomLabel.frame     = NSMakeRect(kMargin + 30, zoomY + 4, 64, 18);
    _zoomInButton.frame  = NSMakeRect(kMargin + 98, zoomY, 26, 26);
    _status.frame        = NSMakeRect(kMargin + 136, zoomY + 5, cw - 136, 16);

    _presetLabel.frame = NSMakeRect(kMargin, zoomY - 40, 50, 17);
    _presets.frame     = NSMakeRect(kMargin + 54, zoomY - 45, std::min<CGFloat>(360, cw - 54), 25);
    _presetDesc.frame  = NSMakeRect(kMargin + 54, zoomY - 66, cw - 54, 16);
}

- (void)windowDidResize:(NSNotification*)n {
    [self layoutViews];
    // Refetching from Photoshop on every frame of a live drag would crawl, so
    // wait until the size settles.
    [NSObject cancelPreviousPerformRequestsWithTarget:self
                                             selector:@selector(previewSizeSettled)
                                               object:nil];
    [self performSelector:@selector(previewSizeSettled) withObject:nil afterDelay:0.15];
}

- (void)previewSizeSettled {
    if (![self refetch]) return;
    [self updatePannable];
    [self rerender:nil];
}

// -------------------------------------------------------------------- values --

// Slider order, written down once:
//   0 squeeze  1 distortion  2 vignette  3 chromatic aberration  4 aperture
//   5 edge softness  6 astigmatism  7 coma  8 grain  9 grain size  10 grain colour
- (void)readControls {
    _controls->squeeze             = std::round(_sliders[0].doubleValue * 100.0) / 100.0;
    _controls->distortion          = std::round(_sliders[1].doubleValue);
    _controls->vignette            = std::round(_sliders[2].doubleValue);
    _controls->chromaticAberration = std::round(_sliders[3].doubleValue);
    _controls->tStop               = std::round(_sliders[4].doubleValue * 10.0) / 10.0;
    _controls->edgeBlur            = std::round(_sliders[5].doubleValue);
    _controls->astigmatism         = std::round(_sliders[6].doubleValue);
    _controls->coma                = std::round(_sliders[7].doubleValue);
    _controls->grain               = std::round(_sliders[8].doubleValue);
    _controls->grainSize           = std::round(_sliders[9].doubleValue * 10.0) / 10.0;
    _controls->grainColour         = std::round(_sliders[10].doubleValue);
    _controls->bands = (_qBest.state == NSControlStateValueOn)  ? 15
                     : (_qDraft.state == NSControlStateValueOn) ?  7 : 11;
}

- (void)syncFields {
    // decimals per row: whole numbers except aperture, squeeze and grain size
    static const int dp[11] = {2, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
    for (int i = 0; i < 11; ++i) {
        if (dp[i] == 0)
            _fields[i].stringValue = [NSString stringWithFormat:@"%d", int(_sliders[i].doubleValue)];
        else
            _fields[i].stringValue =
                [NSString stringWithFormat:@"%.*f", dp[i], _sliders[i].doubleValue];
    }
}

- (void)pushToSliders:(const LensControls&)c {
    _sliders[0].doubleValue = c.squeeze;
    _sliders[1].doubleValue = c.distortion;
    _sliders[2].doubleValue = c.vignette;
    _sliders[3].doubleValue = c.chromaticAberration;
    _sliders[4].doubleValue = c.tStop;
    _sliders[5].doubleValue = c.edgeBlur;
    _sliders[6].doubleValue = c.astigmatism;
    _sliders[7].doubleValue = c.coma;
    _sliders[8].doubleValue = c.grain;
    _sliders[9].doubleValue = c.grainSize;
    _sliders[10].doubleValue = c.grainColour;
    _qDraft.state  = (c.bands <= 7)  ? NSControlStateValueOn : NSControlStateValueOff;
    _qNormal.state = (c.bands > 7 && c.bands < 15) ? NSControlStateValueOn : NSControlStateValueOff;
    _qBest.state   = (c.bands >= 15) ? NSControlStateValueOn : NSControlStateValueOff;
    [self readControls];
    [self syncFields];
    [self rerender:nil];
}

static NSString* const kSaveItem   = @"Save current settings…";
static NSString* const kDeleteItem = @"Delete a saved preset…";
static NSString* const kExportItem = @"Export current settings to a file…";
static NSString* const kImportItem = @"Load a preset from a file…";

- (void)rebuildPresetMenu {
    NSString* keep = _presets.titleOfSelectedItem;
    [_presets removeAllItems];

    std::string lastGroup;
    for (const auto& pr : lensPresets()) {
        if (std::string(pr.group) != lastGroup && !lastGroup.empty())
            [_presets.menu addItem:[NSMenuItem separatorItem]];
        lastGroup = pr.group;
        [_presets addItemWithTitle:@(pr.name)];
    }

    NSDictionary* saved = [presetStore() dictionaryForKey:kPresetsKey];
    NSArray* names = [saved.allKeys sortedArrayUsingSelector:@selector(caseInsensitiveCompare:)];
    if (names.count > 0) {
        [_presets.menu addItem:[NSMenuItem separatorItem]];
        for (NSString* n in names) [_presets addItemWithTitle:n];
    }

    [_presets.menu addItem:[NSMenuItem separatorItem]];
    [_presets addItemWithTitle:kSaveItem];
    if (names.count > 0) [_presets addItemWithTitle:kDeleteItem];
    [_presets addItemWithTitle:kExportItem];
    [_presets addItemWithTitle:kImportItem];

    if (keep && [_presets indexOfItemWithTitle:keep] >= 0) [_presets selectItemWithTitle:keep];
    else [_presets selectItemAtIndex:0];
}

- (void)saveCurrent {
    [self readControls];

    NSAlert* a = [[NSAlert alloc] init];
    a.messageText = @"Save these settings as a preset";
    a.informativeText = @"Saved presets appear in this menu from now on.";
    [a addButtonWithTitle:@"Save"];
    [a addButtonWithTitle:@"Cancel"];
    NSTextField* field = [[NSTextField alloc] initWithFrame:NSMakeRect(0, 0, 240, 24)];
    field.placeholderString = @"My lens";
    a.accessoryView = field;
    [a.window setInitialFirstResponder:field];

    if ([a runModal] != NSAlertFirstButtonReturn) return;
    NSString* name = [field.stringValue stringByTrimmingCharactersInSet:
                          [NSCharacterSet whitespaceAndNewlineCharacterSet]];
    if (name.length == 0) return;

    NSUserDefaults* store = presetStore();
    NSMutableDictionary* saved =
        [([store dictionaryForKey:kPresetsKey] ?: @{}) mutableCopy];
    saved[name] = controlsToDict(*_controls);
    [store setObject:saved forKey:kPresetsKey];
    [store synchronize];

    [self rebuildPresetMenu];
    [_presets selectItemWithTitle:name];
}

- (void)deleteSaved {
    NSUserDefaults* store = presetStore();
    NSMutableDictionary* saved =
        [([store dictionaryForKey:kPresetsKey] ?: @{}) mutableCopy];
    NSArray* names = [saved.allKeys sortedArrayUsingSelector:@selector(caseInsensitiveCompare:)];
    if (names.count == 0) return;

    NSAlert* a = [[NSAlert alloc] init];
    a.messageText = @"Delete a saved preset";
    [a addButtonWithTitle:@"Delete"];
    [a addButtonWithTitle:@"Cancel"];
    NSPopUpButton* pick = [[NSPopUpButton alloc] initWithFrame:NSMakeRect(0, 0, 240, 25)];
    [pick addItemsWithTitles:names];
    a.accessoryView = pick;

    if ([a runModal] != NSAlertFirstButtonReturn) return;
    [saved removeObjectForKey:pick.titleOfSelectedItem];
    [store setObject:saved forKey:kPresetsKey];
    [store synchronize];
    [self rebuildPresetMenu];
}

- (void)presetChosen:(id)sender {
    // Match by TITLE, not by index: the menu carries separators and command items
    // that the preset list does not, so the two indices do not line up.
    NSString* chosen = _presets.titleOfSelectedItem;

    if ([chosen isEqualToString:kSaveItem])   { [self saveCurrent]; return; }
    if ([chosen isEqualToString:kDeleteItem]) { [self deleteSaved]; return; }
    if ([chosen isEqualToString:kExportItem]) { [self exportToFile]; return; }
    if ([chosen isEqualToString:kImportItem]) { [self importFromFile]; return; }

    for (const auto& pr : lensPresets())
        if ([chosen isEqualToString:@(pr.name)]) {
            [self pushToSliders:pr.c];
            self->_presetDesc.stringValue = @(pr.desc);
            return;
        }

    NSDictionary* saved = [presetStore() dictionaryForKey:kPresetsKey];
    NSDictionary* one = saved[chosen];
    if ([one isKindOfClass:[NSDictionary class]]) {
        [self pushToSliders:dictToControls(one)];
        _presetDesc.stringValue = @"Your own saved preset.";
    }
}

// Presets travel as a small JSON file, so they can be shared or kept in a
// project folder rather than living only in this machine's preferences.
- (void)exportToFile {
    [self readControls];
    NSSavePanel* sp = [NSSavePanel savePanel];
    sp.nameFieldStringValue = @"lens.filmic.json";
    sp.allowedFileTypes = @[@"json"];
    if ([sp runModal] != NSModalResponseOK || sp.URL == nil) return;

    NSError* err = nil;
    NSData* data = [NSJSONSerialization dataWithJSONObject:controlsToDict(*_controls)
                                                    options:NSJSONWritingPrettyPrinted
                                                      error:&err];
    if (data == nil || ![data writeToURL:sp.URL atomically:YES]) {
        _presetDesc.stringValue = @"Could not write that file.";
        return;
    }
    _presetDesc.stringValue = [@"Saved to " stringByAppendingString:sp.URL.lastPathComponent];
}

- (void)importFromFile {
    NSOpenPanel* op = [NSOpenPanel openPanel];
    op.allowedFileTypes = @[@"json"];
    op.allowsMultipleSelection = NO;
    if ([op runModal] != NSModalResponseOK || op.URL == nil) return;

    NSData* data = [NSData dataWithContentsOfURL:op.URL];
    id obj = data ? [NSJSONSerialization JSONObjectWithData:data options:0 error:nil] : nil;
    if (![obj isKindOfClass:[NSDictionary class]]) {
        _presetDesc.stringValue = @"That file is not a Filmic preset.";
        return;
    }
    [self pushToSliders:dictToControls(obj)];
    _presetDesc.stringValue = [@"Loaded " stringByAppendingString:op.URL.lastPathComponent];
}

- (void)sliderMoved:(NSSlider*)s {
    // A continuous slider sends its value all through the drag and once more on
    // mouse-up. Render coarse while the mouse is down and properly when it is
    // released: a preview that trails the slider by a second is worse than one
    // that is briefly rough.
    const NSEvent* e = [NSApp currentEvent];
    _dragging = (e != nil && e.type != NSEventTypeLeftMouseUp);
    [self readControls];
    [self syncFields];
    [self rerender:nil];
}

- (void)fieldEdited:(NSTextField*)f {
    NSSlider* s = _sliders[f.tag];
    double v = std::clamp(f.doubleValue, s.minValue, s.maxValue);
    s.doubleValue = v;
    [self readControls];
    [self syncFields];
    [self rerender:nil];
}

- (void)reset:(id)sender {
    [_presets selectItemAtIndex:0];
    [self pushToSliders:LensControls{}];
}

- (void)ok:(id)sender     { [self readControls]; self.accepted = YES; [NSApp stopModal]; }
- (void)cancel:(id)sender { self.accepted = NO;  [NSApp stopModal]; }

// Option is read HERE, at the moment of the click, rather than being baked into
// the button's action ahead of time.
- (void)cancelOrReset:(id)sender {
    if (([NSEvent modifierFlags] & NSEventModifierFlagOption) != 0) [self reset:sender];
    else                                                            [self cancel:sender];
}

// ---------------------------------------------------------------------- zoom --

- (void)updateZoomLabel {
    _zoomLabel.stringValue =
        [NSString stringWithFormat:@"%d%%", int(std::lround(100.0 / double(_zoomScale)))];
}

- (void)zoomTo:(int)scale {
    scale = std::clamp(scale, 1, std::max(1, _fitScale));
    if (scale == _zoomScale) return;
    _zoomScale = scale;
    if (![self refetch]) return;
    [self updateZoomLabel];
    [self updatePannable];
    [self rerender:nil];
}

- (void)zoomIn:(id)sender  { [self zoomTo:_zoomScale / 2]; }
- (void)zoomOut:(id)sender { [self zoomTo:_zoomScale * 2]; }

// How big the whole filter area is at the current zoom, in preview pixels. When
// that fits inside the well there is nothing to pan to.
- (void)fullProxyWidth:(int*)fw height:(int*)fh {
    int aw = 0, ah = 0;
    _source->areaSize(aw, ah);
    *fw = (aw + _zoomScale - 1) / _zoomScale;
    *fh = (ah + _zoomScale - 1) / _zoomScale;
}

- (void)updatePannable {
    int fw = 0, fh = 0;
    [self fullProxyWidth:&fw height:&fh];
    _preview.pannable = (fw > _pw) || (fh > _ph);
    [self.window invalidateCursorRectsForView:_preview];
}

// Dragging moves the picture with the pointer, so the view centre moves the
// other way. One view pixel is one preview pixel: the well and the proxy are
// sized to match.
- (void)panByX:(CGFloat)dx y:(CGFloat)dy {
    int fw = 0, fh = 0;
    [self fullProxyWidth:&fw height:&fh];
    if (fw <= _pw && fh <= _ph) return;

    const double before[2] = {_centreX, _centreY};
    _centreX = std::clamp(_centreX - double(dx) / double(fw), 0.0, 1.0);
    // The view's y runs up, the image's runs down.
    _centreY = std::clamp(_centreY + double(dy) / double(fh), 0.0, 1.0);
    if (_centreX == before[0] && _centreY == before[1]) return;

    if (![self refetch]) return;
    [self rerender:nil];
}

// ------------------------------------------------------------------- preview --

// A grid of straight lines and a lattice of dots, laid into the picture BEFORE
// the optics run. A percentage is abstract; a straight line that bows, fringes
// and softens is not. This is why the reading is worth having.
- (void)drawGuidesInto:(std::vector<float>&)buf {
    const int step = std::max(24, _pw / 12);
    auto put = [&](int x, int y, float v) {
        if (x < 0 || y < 0 || x >= _pw || y >= _ph) return;
        const size_t o = (size_t(y) * _pw + x) * 3;
        buf[o] = buf[o + 1] = buf[o + 2] = v;
    };
    for (int y = 0; y < _ph; ++y)
        for (int x = step / 2; x < _pw; x += step) put(x, y, 0.95f);
    for (int x = 0; x < _pw; ++x)
        for (int y = step / 2; y < _ph; y += step) put(x, y, 0.95f);
    // Dots at every crossing: these show the shape of the softness.
    for (int y = step / 2; y < _ph; y += step)
        for (int x = step / 2; x < _pw; x += step)
            for (int dy = -2; dy <= 2; ++dy)
                for (int dx = -2; dx <= 2; ++dx)
                    if (dx * dx + dy * dy <= 4) put(x + dx, y + dy, 1.0f);
}

- (NSImage*)imageFrom:(const std::vector<float>&)rgb width:(int)_w height:(int)_h {
    NSBitmapImageRep* rep = [[NSBitmapImageRep alloc]
        initWithBitmapDataPlanes:NULL pixelsWide:_w pixelsHigh:_h
                   bitsPerSample:8 samplesPerPixel:3 hasAlpha:NO isPlanar:NO
                  colorSpaceName:NSDeviceRGBColorSpace
                     bytesPerRow:_w * 3 bitsPerPixel:24];
    unsigned char* out = [rep bitmapData];
    const size_t n = size_t(_w) * _h * 3;
    for (size_t i = 0; i < n; ++i) {
        // Already in the document's own encoding, scaled to 0..1 on the way in,
        // so this only undoes that scaling. No extra gamma.
        const float v = rgb[i] * 255.0f;
        out[i] = (unsigned char)(v < 0 ? 0 : (v > 255 ? 255 : v + 0.5f));
    }
    NSImage* img = [[NSImage alloc] initWithSize:NSMakeSize(_w, _h)];
    [img addRepresentation:rep];
    return img;
}

- (void)rerender:(id)sender {
    [self readControls];
    if (_renderRunning) { _renderPending = YES; return; }

    std::vector<float> base = _original;
    if (_gridToggle.state == NSControlStateValueOn) [self drawGuidesInto:base];

    // Dragging renders at the SAME size and fewer wavelengths, not at half
    // resolution. Halving the resolution was four times cheaper but swapped the
    // image for a softer one on every drag and snapped it back on release, which
    // reads as the preview glitching. Band count is the one thing that can be
    // cut without the picture visibly changing shape or sharpness: it only
    // affects how finely colour is resolved in the fringing.
    const int rw = _pw, rh = _ph;
    const double rscale = _pixelScale;

    const BOOL showResult = (_previewToggle.state == NSControlStateValueOn);
    LensControls c = *_controls;
    if (_dragging) c.bands = 7;    // Draft while the mouse is down

    if (!showResult || !c.any()) {
        _preview.image = [self imageFrom:base width:_pw height:_ph];
        _status.stringValue = showResult ? @"Every control is at zero."
                                         : @"Preview off — showing the original.";
        return;
    }

    _renderRunning = YES;
    _status.stringValue = @"Rendering…";

    const int w = rw, h = rh;
    const double scale = rscale;
    const std::string table = _tablePath;
    const BOOL coarse = _dragging;
    __block std::vector<float> buffer = std::move(base);

    dispatch_async(dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^{
        const NSTimeInterval t0 = [NSDate timeIntervalSinceReferenceDate];
        const std::string err = lensRender(buffer.data(), w, h, c, scale, table);
        const int ms = int(([NSDate timeIntervalSinceReferenceDate] - t0) * 1000.0);
        dispatch_async(dispatch_get_main_queue(), ^{
            if (err.empty()) {
                self->_preview.image = [self imageFrom:buffer width:w height:h];
                self->_status.stringValue = coarse
                    ? [NSString stringWithFormat:@"%d × %d draft, %d ms", w, h, ms]
                    : [NSString stringWithFormat:@"%d × %d preview, %d ms", w, h, ms];
            } else {
                self->_status.stringValue =
                    [NSString stringWithUTF8String:("Preview failed: " + err).c_str()];
            }
            self->_renderRunning = NO;
            // A coarse frame finishing after the drag ended must be followed by
            // the proper one, or the preview stays rough until the next nudge.
            if (self->_renderPending || (coarse && !self->_dragging)) {
                self->_renderPending = NO;
                [self rerender:nil];
            }
        });
    });
}

@end

// ---------------------------------------------------------------------------

bool showLensDialog(LensControls& controls, ProxySource& source) {
    @autoreleasepool {
        LensDialogController* c =
            [[LensDialogController alloc] initWithControls:&controls source:&source];
        if (c == nil) {
            // A plug-in that opens nothing and says nothing is the worst possible
            // failure to diagnose.
            NSAlert* a = [[NSAlert alloc] init];
            a.messageText = @"Filmic could not read the image";
            a.informativeText = source.lastError.empty()
                ? @"Photoshop did not supply any pixels to preview."
                : [NSString stringWithUTF8String:source.lastError.c_str()];
            [a runModal];
            return false;
        }
        [NSApp runModalForWindow:c.window];
        [c.window orderOut:nil];
        return c.accepted ? true : false;
    }
}
