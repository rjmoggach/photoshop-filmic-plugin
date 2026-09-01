// Filmic Lens -- the UXP half. It owns the controls and the pixels; all optics
// happen in the native addon.
//
// Flow: read the active layer through the Imaging API, hand the float buffer to
// lensaddon, write the result back inside a history-undoable modal scope.

const { app, core, imaging, action } = require("photoshop");
const fs = require("uxp").storage.localFileSystem;

let addon = null;
let tableToken = null;

const CONTROLS = [
  "chromaticAberration", "vignette", "edgeBlur",
  "astigmatism", "coma", "tStop", "bands",
];

function say(msg, isError) {
  const el = document.getElementById("status");
  el.textContent = msg;
  el.className = isError ? "err" : "";
  if (isError) console.error(msg);
}

function label(id, raw) {
  if (id === "tStop") return "T/" + (raw / 10).toFixed(1);
  if (id === "bands") return String(raw);
  return raw + "%";
}

function readControls() {
  const p = {};
  for (const id of CONTROLS) p[id] = Number(document.getElementById(id).value);
  p.tStop = p.tStop / 10;          // slider is in tenths so it can step finely
  p.tStopWide = 2.0;
  p.highlightRecovery = 1;
  return p;
}

// The spectral table ships beside the plugin. The addon reads it once and caches.
async function tablePath() {
  if (tableToken) return tableToken;
  const folder = await fs.getPluginFolder();
  const entry = await folder.getEntry("data/rec2020.bin");
  tableToken = { path: entry.nativePath };
  return tableToken;
}

async function applyToLayer() {
  const doc = app.activeDocument;
  if (!doc) { say("Open a document first.", true); return; }
  const layer = doc.activeLayers && doc.activeLayers[0];
  if (!layer) { say("Select a layer first.", true); return; }

  const params = readControls();
  const anyEffect = params.chromaticAberration || params.vignette ||
                    params.edgeBlur || params.astigmatism || params.coma;
  if (!anyEffect) { say("Every control is at zero — nothing to apply.", true); return; }

  say("Reading pixels…");
  const t0 = Date.now();

  await core.executeAsModal(async (ctx) => {
    ctx.reportProgress({ value: 0.05 });

    const px = await imaging.getPixels({
      documentID: doc.id,
      layerID: layer.id,
      componentSize: 32,          // float32 per channel
      applyAlpha: false,
    });

    const w = px.imageData.width;
    const h = px.imageData.height;
    const comps = px.imageData.components;
    const raw = await px.imageData.getData({ chunky: true });

    // lenscore works in interleaved RGB. Drop alpha on the way in and restore it
    // on the way out, so a layer mask or transparency survives untouched.
    const rgb = new Float32Array(w * h * 3);
    for (let i = 0, s = 0, d = 0; i < w * h; i++, s += comps, d += 3) {
      rgb[d] = raw[s]; rgb[d + 1] = raw[s + 1]; rgb[d + 2] = raw[s + 2];
    }

    ctx.reportProgress({ value: 0.2 });
    say(`Rendering ${w}×${h}…`);

    const table = await tablePath();
    // The addon reads the ArrayBuffer directly -- this SDK has no
    // typed-array accessor, only get_arraybuffer_info.
    const out = addon.renderPixels(rgb.buffer, w, h, params, table);

    ctx.reportProgress({ value: 0.8 });

    const result = new Float32Array(raw.length);
    for (let i = 0, s = 0, d = 0; i < w * h; i++, s += 3, d += comps) {
      result[d] = out[s]; result[d + 1] = out[s + 1]; result[d + 2] = out[s + 2];
      for (let c = 3; c < comps; c++) result[d + c] = raw[d + c];   // alpha through
    }

    const outData = await imaging.createImageDataFromBuffer(result, {
      width: w, height: h, components: comps,
      colorSpace: px.imageData.colorSpace, componentSize: 32, chunky: true,
    });

    await imaging.putPixels({ documentID: doc.id, layerID: layer.id, imageData: outData });

    px.imageData.dispose();
    outData.dispose();
    ctx.reportProgress({ value: 1 });
  }, { commandName: "Filmic Lens", interactive: false });

  say(`Done in ${((Date.now() - t0) / 1000).toFixed(1)}s.`);
}

function wireUp() {
  for (const id of CONTROLS) {
    const el = document.getElementById(id);
    const out = document.getElementById(id + "-val");
    out.textContent = label(id, Number(el.value));
    el.addEventListener("input", () => { out.textContent = label(id, Number(el.value)); });
  }

  document.getElementById("apply").addEventListener("click", () => {
    applyToLayer().catch((e) => say("Failed: " + (e && e.message ? e.message : e), true));
  });

  try {
    addon = require("lens.uxpaddon");
    say(addon.version ? addon.version() + " ready." : "Addon loaded.");
  } catch (e) {
    say("Native addon failed to load — is lens.uxpaddon built and staged? " +
        (e && e.message ? e.message : e), true);
  }
}

document.addEventListener("DOMContentLoaded", wireUp);
if (document.readyState !== "loading") wireUp();
