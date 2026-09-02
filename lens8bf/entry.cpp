// Filmic -- a Photoshop filter plug-in (8BF), under Filter > Moggach > Filmic.
//
// Unlike the UXP panel this replaces, a filter runs inside Photoshop's own
// filter machinery: it gets a live proxy of the image for its dialog, it writes
// back over the pixels it was given rather than creating a layer, and applied to
// a Smart Object it becomes a re-editable Smart Filter for free.
//
// The whole selection is filtered in one pass. That is not laziness about
// tiling: the optics are field-dependent -- every kernel depends on where in the
// FRAME its pixel sits -- so a tile cannot be rendered correctly without knowing
// the whole frame anyway.

#include "PIDefines.h"
#include "PITypes.h"
#include "PIFilter.h"
#include "PIAbout.h"
#include "PIGeneral.h"

#include "LensParams.h"
#include "dialog.h"
#include "render.h"

#include <algorithm>
#include <string>
#include <vector>

#include <CoreFoundation/CoreFoundation.h>

FilterRecord* gFilterRecord = nullptr;
SPBasicSuite* sSPBasic      = nullptr;
static int16* gResult       = nullptr;   // Photoshop passes this alongside the record

namespace {

LensControls gControls;

// ---------------------------------------------------------------- rectangles --
// Photoshop keeps two sets of rectangle fields: the historical 16-bit ones and
// the 32-bit "big document" ones. Which set is live depends on whether the host
// handed us bigDocumentData, so every read and write goes through these.

void setInRect(const VRect& r) {
    if (gFilterRecord->bigDocumentData != nullptr) {
        gFilterRecord->bigDocumentData->inRect32 = r;
        gFilterRecord->bigDocumentData->PluginUsing32BitCoordinates = true;
    } else {
        gFilterRecord->inRect.top = int16(r.top);       gFilterRecord->inRect.left  = int16(r.left);
        gFilterRecord->inRect.bottom = int16(r.bottom); gFilterRecord->inRect.right = int16(r.right);
    }
}

void setOutRect(const VRect& r) {
    if (gFilterRecord->bigDocumentData != nullptr) {
        gFilterRecord->bigDocumentData->outRect32 = r;
        gFilterRecord->bigDocumentData->PluginUsing32BitCoordinates = true;
    } else {
        gFilterRecord->outRect.top = int16(r.top);       gFilterRecord->outRect.left  = int16(r.left);
        gFilterRecord->outRect.bottom = int16(r.bottom); gFilterRecord->outRect.right = int16(r.right);
    }
}

void setMaskRect(const VRect& r) {
    if (gFilterRecord->bigDocumentData != nullptr) {
        gFilterRecord->bigDocumentData->maskRect32 = r;
        gFilterRecord->bigDocumentData->PluginUsing32BitCoordinates = true;
    } else {
        gFilterRecord->maskRect.top = int16(r.top);       gFilterRecord->maskRect.left  = int16(r.left);
        gFilterRecord->maskRect.bottom = int16(r.bottom); gFilterRecord->maskRect.right = int16(r.right);
    }
}

VRect filterRect() {
    VRect r{};
    if (gFilterRecord->bigDocumentData != nullptr) {
        r = gFilterRecord->bigDocumentData->filterRect32;
    } else {
        r.top = gFilterRecord->filterRect.top;       r.left  = gFilterRecord->filterRect.left;
        r.bottom = gFilterRecord->filterRect.bottom; r.right = gFilterRecord->filterRect.right;
    }
    return r;
}

const VRect kEmptyRect{0, 0, 0, 0};

// -------------------------------------------------------------------- pixels --
// Photoshop hands over 8, 16 or 32 bits per channel. lenscore works in float
// 0..1, so one scale converts in each direction -- the SAME constant both ways,
// which keeps the round trip exact whatever Photoshop's white point for a depth
// turns out to be.

double depthScale(int16 depth) {
    if (depth == 8)  return 255.0;
    if (depth == 16) return 32768.0;   // Photoshop's 16-bit white
    return 1.0;                        // 32-bit is already float 0..1
}

int planeCount() {
    // Filter colour only. Alpha and any extra channels are left untouched.
    return std::min<int>(gFilterRecord->planes, gFilterRecord->imageMode == plugInModeGrayScale ||
                                                gFilterRecord->imageMode == plugInModeGray16 ? 1 : 3);
}

// Copy Photoshop's buffer into interleaved float RGB. Grayscale is widened to
// three equal channels so one code path renders everything.
void readPixels(const void* data, int32 rowBytes, int w, int h, int planes,
                int16 depth, std::vector<float>& rgb) {
    rgb.assign(size_t(w) * h * 3, 0.0f);
    const double inv = 1.0 / depthScale(depth);
    const auto* base = static_cast<const uint8*>(data);

    for (int y = 0; y < h; ++y) {
        const uint8* row = base + size_t(y) * rowBytes;
        for (int x = 0; x < w; ++x) {
            float v[3] = {0, 0, 0};
            for (int c = 0; c < planes; ++c) {
                double s = 0.0;
                if (depth == 8)       s = reinterpret_cast<const uint8*>(row)[size_t(x) * planes + c];
                else if (depth == 16) s = reinterpret_cast<const uint16*>(row)[size_t(x) * planes + c];
                else                  s = reinterpret_cast<const float*>(row)[size_t(x) * planes + c];
                v[c] = float(s * inv);
            }
            if (planes == 1) { v[1] = v[0]; v[2] = v[0]; }
            const size_t o = (size_t(y) * w + x) * 3;
            rgb[o] = v[0]; rgb[o + 1] = v[1]; rgb[o + 2] = v[2];
        }
    }
}

void writePixels(void* data, int32 rowBytes, int w, int h, int planes,
                 int16 depth, const std::vector<float>& rgb) {
    const double scale = depthScale(depth);
    const double top   = (depth == 32) ? 0.0 : scale;   // 0 means "do not clamp"
    auto* base = static_cast<uint8*>(data);

    for (int y = 0; y < h; ++y) {
        uint8* row = base + size_t(y) * rowBytes;
        for (int x = 0; x < w; ++x) {
            const size_t o = (size_t(y) * w + x) * 3;
            for (int c = 0; c < planes; ++c) {
                double v = double(rgb[o + (planes == 1 ? 0 : c)]) * scale;
                if (top > 0.0) v = std::clamp(v, 0.0, top);
                if (depth == 8)       reinterpret_cast<uint8*>(row)[size_t(x) * planes + c]  = uint8(v + 0.5);
                else if (depth == 16) reinterpret_cast<uint16*>(row)[size_t(x) * planes + c] = uint16(v + 0.5);
                else                  reinterpret_cast<float*>(row)[size_t(x) * planes + c]  = float(v);
            }
        }
    }
}

// ------------------------------------------------------------- table on disk --
// The spectral table sits beside the executable inside the plug-in bundle.
std::string tablePathInBundle() {
    std::string path;
    CFBundleRef me = CFBundleGetBundleWithIdentifier(CFSTR("com.moggach.filmic"));
    if (me == nullptr) return path;
    CFURLRef url = CFBundleCopyResourceURL(me, CFSTR("rec2020"), CFSTR("bin"), nullptr);
    if (url == nullptr) return path;
    char buf[PATH_MAX] = {0};
    if (CFURLGetFileSystemRepresentation(url, true, reinterpret_cast<UInt8*>(buf), PATH_MAX))
        path = buf;
    CFRelease(url);
    return path;
}

// ------------------------------------------------------------- proxy source --
// How the dialog gets something to draw. Photoshop will scale the image down for
// us: set inputRate and ask for the whole filter rectangle, and advanceState
// hands back a proxy at that scale.
class PhotoshopProxy : public ProxySource {
public:
    bool fetchProxy(int maxW, int maxH, int scale, double cxFrac, double cyFrac,
                    std::vector<float>& rgb, int& w, int& h, double& pixelScale) override {
        const VRect r = filterRect();
        const int32 fw = r.right - r.left, fh = r.bottom - r.top;
        if (fw <= 0 || fh <= 0) { lastError = "the filter area is empty"; return false; }

        if (scale < 1) {                       // 0 asks for "fit the whole area"
            // Fit BOTH dimensions into the well, which is not square.
            scale = int(std::max((fw + maxW - 1) / maxW, (fh + maxH - 1) / maxH));
            if (scale < 1) scale = 1;
        }

        // The whole area, expressed in proxy pixels, and the window of it we
        // actually want. inRect is given in PROXY coordinates, not
        // full-resolution ones -- the rectangle is scaled down by the same
        // factor as inputRate. Full-resolution coordinates here ask for an area
        // far outside the proxy image and advanceState refuses the request.
        const int32 pw = (fw + scale - 1) / scale;
        const int32 ph = (fh + scale - 1) / scale;
        const int32 vw = std::min<int32>(maxW, pw);
        const int32 vh = std::min<int32>(maxH, ph);

        int32 left = int32(cxFrac * double(pw) - double(vw) * 0.5);
        int32 top  = int32(cyFrac * double(ph) - double(vh) * 0.5);
        left = std::clamp<int32>(left, 0, pw - vw);
        top  = std::clamp<int32>(top,  0, ph - vh);

        VRect proxy{};
        proxy.left   = r.left / scale + left;
        proxy.top    = r.top  / scale + top;
        proxy.right  = proxy.left + vw;
        proxy.bottom = proxy.top  + vh;
        if (vw <= 0 || vh <= 0) { lastError = "the filter area is too small to preview"; return false; }

        setInRect(proxy);
        setMaskRect(proxy);
        setOutRect(kEmptyRect);
        gFilterRecord->inputRate    = int32(scale) << 16;
        gFilterRecord->maskRate     = int32(scale) << 16;
        gFilterRecord->inputPadding = 255;
        gFilterRecord->maskPadding  = 255;
        gFilterRecord->inLoPlane    = 0;
        gFilterRecord->inHiPlane    = int16(planeCount() - 1);

        const OSErr err = gFilterRecord->advanceState();
        if (err != noErr) {
            lastError = "Photoshop refused the preview request (error " +
                        std::to_string(int(err)) + ")";
            return false;
        }
        if (gFilterRecord->inData == nullptr) {
            lastError = "Photoshop returned no preview pixels";
            return false;
        }

        w = int(vw);
        h = int(vh);
        readPixels(gFilterRecord->inData, gFilterRecord->inRowBytes, w, h,
                   planeCount(), gFilterRecord->depth, rgb);
        pixelScale = double(scale);
        return true;
    }

    void areaSize(int& w, int& h) const override {
        const VRect r = filterRect();
        w = int(r.right - r.left);
        h = int(r.bottom - r.top);
    }

    std::string tablePath() const override { return tablePathInBundle(); }
};

// ----------------------------------------------------------------- selectors --

void doStart() {
    PhotoshopProxy proxy;

    // The dialog owns the preview loop; it calls back into the proxy above.
    if (!showLensDialog(gControls, proxy)) {
        *gResult = userCanceledErr;
        return;
    }
    if (!gControls.any()) return;   // nothing asked for, nothing to do

    const VRect r = filterRect();
    const int w = int(r.right - r.left), h = int(r.bottom - r.top);
    if (w <= 0 || h <= 0) return;

    // Full resolution this time.
    gFilterRecord->inputRate  = int32(1) << 16;
    gFilterRecord->maskRate   = int32(1) << 16;
    gFilterRecord->inLoPlane  = 0;
    gFilterRecord->inHiPlane  = int16(planeCount() - 1);
    gFilterRecord->outLoPlane = 0;
    gFilterRecord->outHiPlane = int16(planeCount() - 1);
    setInRect(r);
    setOutRect(r);
    setMaskRect(kEmptyRect);

    if (gFilterRecord->advanceState() != noErr) { *gResult = filterBadParameters; return; }
    if (gFilterRecord->inData == nullptr || gFilterRecord->outData == nullptr) return;

    std::vector<float> rgb;
    readPixels(gFilterRecord->inData, gFilterRecord->inRowBytes, w, h,
               planeCount(), gFilterRecord->depth, rgb);

    const std::string err = lensRender(rgb.data(), w, h, gControls, 1.0, tablePathInBundle());
    if (!err.empty()) { *gResult = filterBadParameters; return; }

    writePixels(gFilterRecord->outData, gFilterRecord->outRowBytes, w, h,
                planeCount(), gFilterRecord->depth, rgb);

    // Nothing left to ask for; Photoshop commits the last outData at Finish.
    setInRect(kEmptyRect);
    setOutRect(kEmptyRect);
    setMaskRect(kEmptyRect);
}

}  // namespace

DLLExport MACPASCAL void PluginMain(const int16 selector, void* filterRecord,
                                    intptr_t* /*data*/, int16* result) {
    if (selector == filterSelectorAbout) {
        sSPBasic = static_cast<AboutRecord*>(filterRecord)->sSPBasic;
        *result = noErr;
        return;
    }

    gFilterRecord = static_cast<FilterRecordPtr>(filterRecord);
    sSPBasic      = gFilterRecord->sSPBasic;
    gResult = result;
    *result = noErr;

    switch (selector) {
        case filterSelectorParameters:
            break;                       // defaults live in LensControls
        case filterSelectorPrepare:
            // The render holds one float RGB copy of the image; ask for nothing
            // from Photoshop's own tile buffers.
            gFilterRecord->bufferSpace = 0;
            gFilterRecord->maxSpace    = 0;
            break;
        case filterSelectorStart:
            doStart();
            break;
        case filterSelectorContinue:
            setInRect(kEmptyRect);
            setOutRect(kEmptyRect);
            setMaskRect(kEmptyRect);
            break;
        case filterSelectorFinish:
            break;
        default:
            *result = filterBadParameters;
            break;
    }
}
