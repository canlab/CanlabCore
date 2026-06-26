import { Niivue, SLICE_TYPE, SHOW_RENDER, NVImage } from './niivue.js';

// canlab_niivue_viewer.js
// ----------------------------------------------------------------------------
// CANlab NiiVue viewer module.
//
// Public API: canlabNiivue(target, config) -> Promise<Niivue>
//
// This module boots a NiiVue volume viewer for CANlab brain-image displays.
// It supports loading an anatomical underlay and a statistical/mask overlay,
// either from URLs or from inline base64-encoded .nii.gz payloads (the form
// the MATLAB writer uses for fully standalone, single-file HTML output).
//
// The module is intentionally dependency-free except for ./niivue.js. The
// MATLAB writer strips the first import line when inlining this module into a
// standalone page (where Niivue/SLICE_TYPE/NVImage are already in scope), so
// that import MUST remain the verbatim first line of the file.
// ----------------------------------------------------------------------------

/** Module version string. */
export const CANLAB_NIIVUE_VERSION = "0.1.0";

/** DOM element ids the viewer wires to (must match the HTML template). */
const CONTROLS_ID = "canlab-niivue-controls";
const READOUT_ID = "canlab-niivue-readout";
/** Opacity-slider panel; created by the viewer if not already in the page. */
const OPACITY_ID = "canlab-niivue-opacity";

/**
 * NiiVue volume colormapType for "colorbar from 0, subthreshold transparent".
 * (COLORMAP_TYPE.ZERO_TO_MAX_TRANSPARENT_BELOW_MIN === 1 in NiiVue 0.57; the enum
 * is not exported from the bundle, so we use the stable numeric literal.)
 * With this set and cal_min > 0, voxels with |value| < cal_min render fully
 * transparent so the anatomical underlay shows through where there is no blob.
 */
const COLORMAP_TYPE_TRANSPARENT_BELOW_MIN = 1;

/**
 * Fixed high-contrast color for the highlighted atlas region (RGB, 0-255).
 * Because only ONE region is ever shown at a time, the whole label LUT is
 * painted this single color rather than a per-region hue, so the highlight
 * reads consistently against grayscale anatomy and hot/cool stat overlays.
 * Bright magenta sits outside both the inferno (warm) and winter (cool) ramps.
 */
const ATLAS_HIGHLIGHT_RGB = [255, 0, 230];

/**
 * Atlas outline thickness, in voxels. NiiVue's label shader marks a voxel as a
 * boundary if a neighbor at +/- this many voxels has a different label, so a
 * larger value yields a thicker outline band. Driven into nv.opts.atlasOutlineWidth
 * (a small CANlab patch in the vendored niivue.js reads it; if that patch is ever
 * lost on a niivue re-download, the outline simply falls back to ~1 voxel).
 */
const ATLAS_OUTLINE_WIDTH = 3;

/** Curated overlay colormaps offered in the dropdown (all exist in NiiVue 0.57). */
const OVERLAY_COLORMAPS = [
    "warm", "hot", "redyellow", "inferno", "plasma", "viridis",
    "turbo", "jet", "cool", "winter", "blue", "gray"
];

/**
 * Map a friendly slice-type name to a NiiVue SLICE_TYPE enum value.
 * @param {string} name - one of multiplanar|axial|coronal|sagittal|render
 * @returns {number} SLICE_TYPE enum value (defaults to MULTIPLANAR)
 */
function sliceTypeFromName(name) {
    switch (String(name || "").toLowerCase()) {
        case "axial":       return SLICE_TYPE.AXIAL;
        case "coronal":     return SLICE_TYPE.CORONAL;
        case "sagittal":    return SLICE_TYPE.SAGITTAL;
        case "render":      return SLICE_TYPE.RENDER;
        case "multiplanar":
        default:            return SLICE_TYPE.MULTIPLANAR;
    }
}

/**
 * Decide whether an image-source config value is a base64 payload or a URL.
 * An object carrying a `.base64` field is treated as a base64 payload; a plain
 * string is treated as a URL.
 * @param {(string|{base64:string,name?:string})} src
 * @returns {boolean} true if `src` is a base64 payload object
 */
function isBase64Source(src) {
    return !!src && typeof src === "object" && typeof src.base64 === "string";
}

/**
 * Load a single image (underlay or overlay) into the NiiVue instance. Both
 * paths are append-style so the underlay (loaded first) is never reloaded and
 * keeps its own options: URL sources go through nv.addVolumeFromUrl, base64
 * sources are decoded via NVImage.loadFromBase64 and pushed with nv.addVolume.
 *
 * @param {Niivue} nv - the NiiVue instance
 * @param {(string|{base64:string,name?:string})} src - image source
 * @param {object} opts - NiiVue per-volume options (colormap, cal_min, ...)
 * @param {string} defaultName - fallback name (must end in .nii.gz)
 * @returns {Promise<void>}
 */
async function loadImage(nv, src, opts, defaultName) {
    if (isBase64Source(src)) {
        // NiiVue detects gzip from a '.nii.gz' name, so require that exact
        // suffix; otherwise fall back to the (.nii.gz) default name.
        const name = src.name && /\.nii\.gz$/i.test(src.name) ? src.name : defaultName;
        const img = await NVImage.loadFromBase64({
            base64: src.base64,
            name,
            ...opts
        });
        nv.addVolume(img);
    } else if (typeof src === "string" && src.length > 0) {
        await nv.addVolumeFromUrl({ url: src, ...opts });
    }
}

/**
 * Build the slice-layout control button bar inside #canlab-niivue-controls.
 * Buttons switch the NiiVue slice type and maintain an `.active` class.
 *
 * @param {Niivue} nv - the NiiVue instance
 * @param {string} initialName - the initial slice-type name to mark active
 * @returns {void}
 */
function buildControls(nv, initialName) {
    const bar = document.getElementById(CONTROLS_ID);
    if (!bar) return;
    bar.innerHTML = "";

    const buttons = [
        { label: "Axial",       name: "axial" },
        { label: "Coronal",     name: "coronal" },
        { label: "Sagittal",    name: "sagittal" },
        { label: "Multiplanar", name: "multiplanar" },
        { label: "3D Render",   name: "render" }
    ];

    const els = [];
    const setActive = (activeName) => {
        for (const el of els) {
            el.classList.toggle("active", el.dataset.slice === activeName);
        }
    };

    for (const spec of buttons) {
        const btn = document.createElement("button");
        btn.type = "button";
        btn.className = "canlab-niivue-btn";
        btn.textContent = spec.label;
        btn.dataset.slice = spec.name;
        btn.addEventListener("click", () => {
            nv.setSliceType(sliceTypeFromName(spec.name));
            setActive(spec.name);
        });
        bar.appendChild(btn);
        els.push(btn);
    }

    setActive(String(initialName || "multiplanar").toLowerCase());
}

/**
 * Build the slider panel (#canlab-niivue-opacity): a dynamic display-threshold
 * slider for the overlay plus overlay/underlay opacity sliders. The panel is
 * created and inserted after the controls bar if not already present, so it
 * works for both the MATLAB-generated template and hand-authored pages.
 *
 * @param {Niivue} nv - the NiiVue instance
 * @param {{hasUnderlay:boolean,hasOverlay:boolean,overlayIndex:number,
 *          overlayOpacity:number,calMin:number,calMax:number}} info
 * @returns {void}
 */
function buildOpacityPanel(nv, info) {
    let panel = document.getElementById(OPACITY_ID);
    if (!panel) {
        panel = document.createElement("div");
        panel.id = OPACITY_ID;
        const controls = document.getElementById(CONTROLS_ID);
        if (controls && controls.parentNode) {
            controls.parentNode.insertBefore(panel, controls.nextSibling);
        } else {
            document.body.insertBefore(panel, document.body.firstChild);
        }
    }
    panel.innerHTML = "";

    // Generalized slider: opts = {min, max, step, value, fmt}.
    const makeSlider = (labelText, opts, onChange) => {
        const fmt = opts.fmt || ((v) => v.toFixed(2));
        const wrap = document.createElement("label");
        wrap.className = "canlab-niivue-slider";

        const name = document.createElement("span");
        name.className = "canlab-niivue-slider-label";
        name.textContent = labelText;

        const input = document.createElement("input");
        input.type = "range";
        input.min = String(opts.min);
        input.max = String(opts.max);
        input.step = String(opts.step);
        input.value = String(opts.value);

        const val = document.createElement("span");
        val.className = "canlab-niivue-slider-value";
        val.textContent = fmt(opts.value);

        input.addEventListener("input", () => {
            const v = parseFloat(input.value);
            val.textContent = fmt(v);
            onChange(v);
        });

        wrap.appendChild(name);
        wrap.appendChild(input);
        wrap.appendChild(val);
        panel.appendChild(wrap);
    };

    // Dynamic display threshold: raise/lower cal_min so |value| < threshold is
    // transparent. Mirrors to the negative colormap so both lobes re-threshold.
    if (info.hasOverlay && info.overlayIndex >= 0 && Number.isFinite(info.calMax)) {
        const tmax = info.calMax;
        const tstart = Number.isFinite(info.calMin) ? info.calMin : 0;
        makeSlider(
            "Threshold",
            { min: 0, max: tmax, step: tmax / 100, value: tstart, fmt: (v) => v.toFixed(2) },
            (v) => {
                const ov = nv.volumes[info.overlayIndex];
                if (!ov) return;
                ov.cal_min = v;
                if ("cal_minNeg" in ov) ov.cal_minNeg = -v;
                nv.updateGLVolume();
            }
        );
    }

    // Overlay then underlay opacity.
    if (info.hasOverlay && info.overlayIndex >= 0) {
        makeSlider("Overlay", { min: 0, max: 1, step: 0.05, value: info.overlayOpacity },
            (v) => nv.setOpacity(info.overlayIndex, v));
    }
    if (info.hasUnderlay) {
        makeSlider("Underlay", { min: 0, max: 1, step: 0.05, value: 1 },
            (v) => nv.setOpacity(0, v));
    }
}

/**
 * Build the "Atlas region" control: a dropdown (Outline / Shaded / Off) plus a
 * controller that highlights ONLY the region under the crosshair. The atlas
 * label volume stays loaded; we show a single region at a time by zeroing every
 * label's alpha in the colormapLabel LUT except the current code, then
 * re-uploading via updateGLVolume(). The LUT mutation is O(nLabels) and the GPU
 * refresh runs only when the crosshair crosses into a different region.
 *
 * @param {Niivue} nv - the NiiVue instance
 * @param {number} atlasIndex - volume index of the atlas label layer
 * @returns {{onCode:(code:number)=>void}|null} controller, or null if unavailable
 */
function buildAtlasControls(nv, atlasIndex) {
    const bar = document.getElementById(CONTROLS_ID);
    const av = nv.volumes && nv.volumes[atlasIndex];
    if (!bar || !av || !av.colormapLabel) return null;

    const lut = av.colormapLabel.lut;            // Uint8ClampedArray, RGBA per dense label
    const lmin = av.colormapLabel.min;
    const lmax = av.colormapLabel.max;
    const nLabel = lmax - lmin + 1;

    // Per-mode look. In NiiVue's label shader the region FILL is drawn at
    // (labelAlpha * opacity), while boundary voxels are drawn at
    // (labelAlpha * atlasOutline), overriding opacity. So an outline-only look
    // needs opacity 0 (no fill) + a high atlasOutline (bright border); a shaded
    // look needs a translucent opacity + atlasOutline 0 (no border emphasis).
    const MODES = {
        outline: { outline: 1.0, opacity: 0.0 },
        shaded:  { outline: 0.0, opacity: 0.5 },
        off:     { outline: 0.0, opacity: 0.0 }
    };
    let mode = "outline";   // default: ON, show the current region's outline
    let curCode = -1;

    // Set LUT alpha so only the current region (curCode) is visible.
    const maskLutToCurrent = () => {
        for (let i = 0; i < nLabel; i++) lut[i * 4 + 3] = 0;
        if (curCode >= 1 && curCode >= lmin && curCode <= lmax) {
            lut[(curCode - lmin) * 4 + 3] = 255;
        }
    };

    const apply = () => {
        const m = MODES[mode] || MODES.off;
        maskLutToCurrent();
        nv.setAtlasOutline(m.outline);          // global, but only the label volume responds
        nv.setOpacity(atlasIndex, mode === "off" ? 0 : m.opacity);
        nv.updateGLVolume();
    };

    // --- Dropdown in the controls bar, right of the "Load overlay" button --
    const wrap = document.createElement("label");
    wrap.className = "canlab-niivue-select";
    wrap.textContent = "Atlas region";

    const sel = document.createElement("select");
    for (const [value, text] of [["outline", "Outline"], ["shaded", "Shaded"], ["off", "Off"]]) {
        const o = document.createElement("option");
        o.value = value;
        o.textContent = text;
        if (value === mode) o.selected = true;
        sel.appendChild(o);
    }
    sel.addEventListener("change", () => {
        mode = sel.value;
        apply();
    });
    wrap.appendChild(sel);
    bar.appendChild(wrap);

    apply();   // establish initial mode (nothing highlighted until a code arrives)

    return {
        onCode(code) {
            if (code === curCode || mode === "off") {
                curCode = code;
                return;
            }
            curCode = code;
            apply();
        }
    };
}

/**
 * Build a NiiVue label colormap (for setColormapLabel) from a region-name list.
 * labels[k] is the name for integer code k+1 (atlas codes are 1..N; code 0 is
 * background, rendered transparent). Every region is painted the single fixed
 * high-contrast highlight color (ATLAS_HIGHLIGHT_RGB) — only one region is shown
 * at a time, so a per-region hue would add no information and reads less clearly.
 *
 * @param {string[]} labels - region names, index = (code - 1)
 * @returns {{R:number[],G:number[],B:number[],A:number[],I:number[],labels:string[]}}
 */
function buildAtlasLabelColormap(labels) {
    const n = Array.isArray(labels) ? labels.length : 0;
    const [hr, hg, hb] = ATLAS_HIGHLIGHT_RGB;
    const R = [0], G = [0], B = [0], A = [0], I = [0], L = [""];
    for (let c = 1; c <= n; c++) {
        R.push(hr); G.push(hg); B.push(hb); A.push(255); I.push(c);
        L.push(String(labels[c - 1] || ""));
    }
    return { R, G, B, A, I, labels: L };
}

/**
 * The current overlay volume is always the last-loaded volume (underlay is index 0).
 * Resolving it dynamically keeps the colormap dropdown correct after a file swap.
 * @param {Niivue} nv
 * @returns {object|null}
 */
function currentOverlay(nv) {
    if (!nv.volumes || nv.volumes.length === 0) return null;
    return nv.volumes[nv.volumes.length - 1];
}

/**
 * Add overlay tools to the controls bar: a colormap dropdown and (optionally) a
 * "Load overlay" file picker that swaps in a user-selected .nii/.nii.gz file.
 *
 * @param {Niivue} nv - the NiiVue instance
 * @param {object} cfg - viewer config (colormap, colormapNegative, opacity, allowLoadOverlay)
 * @returns {void}
 */
function buildOverlayControls(nv, cfg) {
    const bar = document.getElementById(CONTROLS_ID);
    if (!bar) return;

    // --- Colormap dropdown -------------------------------------------------
    const cmapLabel = document.createElement("label");
    cmapLabel.className = "canlab-niivue-select";
    cmapLabel.textContent = "Colormap";

    const sel = document.createElement("select");
    const current = cfg.colormap || "inferno";
    const names = OVERLAY_COLORMAPS.includes(current)
        ? OVERLAY_COLORMAPS
        : [current, ...OVERLAY_COLORMAPS];
    for (const name of names) {
        const opt = document.createElement("option");
        opt.value = name;
        opt.textContent = name;
        if (name === current) opt.selected = true;
        sel.appendChild(opt);
    }
    sel.addEventListener("change", () => {
        const ov = currentOverlay(nv);
        if (ov) {
            nv.setColormap(ov.id, sel.value);
            nv.updateGLVolume();
        }
    });
    cmapLabel.appendChild(sel);
    bar.appendChild(cmapLabel);

    // --- Load-overlay file picker (optional) -------------------------------
    if (cfg.allowLoadOverlay === false) return;

    const fileInput = document.createElement("input");
    fileInput.type = "file";
    fileInput.accept = ".nii,.nii.gz";
    fileInput.style.display = "none";

    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "canlab-niivue-btn";
    btn.textContent = "Load overlay…";
    btn.addEventListener("click", () => fileInput.click());

    fileInput.addEventListener("change", async () => {
        const f = fileInput.files && fileInput.files[0];
        if (!f) return;
        try {
            const img = await NVImage.loadFromFile({ file: f, name: f.name });
            const prev = currentOverlay(nv);
            // Replace the previous overlay (keep the underlay at index 0).
            if (prev && nv.volumes.length > 1) nv.removeVolume(prev);
            nv.addVolume(img);

            const ov = currentOverlay(nv);
            if (ov) {
                nv.setColormap(ov.id, sel.value);
                nv.setColormapNegative(ov.id, cfg.colormapNegative || "winter");
                ov.colormapType = COLORMAP_TYPE_TRANSPARENT_BELOW_MIN;
                nv.setOpacity(nv.volumes.length - 1,
                    cfg.opacity !== undefined && cfg.opacity !== null ? cfg.opacity : 0.8);
                nv.updateGLVolume();
            }
        } catch (e) {
            // eslint-disable-next-line no-console
            console.error("canlab_niivue: failed to load overlay file", e);
        }
        fileInput.value = "";
    });

    bar.appendChild(btn);
    bar.appendChild(fileInput);
}

/**
 * Format a numeric value for the readout, trimming noise but keeping precision.
 * @param {number} v
 * @returns {string}
 */
function fmtValue(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return "n/a";
    if (v === 0) return "0";
    const a = Math.abs(v);
    if (a >= 1000 || a < 0.001) return v.toExponential(3);
    return (Math.round(v * 1000) / 1000).toString();
}

/**
 * Wire the crosshair-location readout into #canlab-niivue-readout. Shows MNI
 * (world) mm coordinates and the value of the overlay at the crosshair, falling
 * back to the underlay value when no overlay is present.
 *
 * @param {Niivue} nv - the NiiVue instance
 * @param {boolean} hasOverlay - whether an overlay volume was loaded
 * @param {{index:number,labels:string[]}} [atlasInfo]
 *        Atlas layer index and label key (labels[k] names integer code k+1).
 *        When present and the crosshair sits on a labeled voxel, the region
 *        name is appended to the readout.
 * @param {{onCode:(code:number)=>void}} [atlasCtl]
 *        Optional atlas-region controller; notified of the integer code under
 *        the crosshair so it can highlight that single region.
 * @returns {void}
 */
function wireReadout(nv, hasOverlay, atlasInfo, atlasCtl) {
    const out = document.getElementById(READOUT_ID);
    if (!out) return;

    out.textContent = "MNI [--- --- ---] mm   value: ---";

    nv.onLocationChange = (data) => {
        if (!data) return;

        // data.mm is a NiiVue vec4 (a Float32Array of [x, y, z, 1]), NOT a plain
        // Array — so test .length, don't use Array.isArray (which is false for typed arrays).
        const mm = data.mm && data.mm.length >= 3 ? data.mm : null;
        const x = mm ? Math.round(mm[0]) : "--";
        const y = mm ? Math.round(mm[1]) : "--";
        const z = mm ? Math.round(mm[2]) : "--";

        // data.values is an array of per-layer {value,...}. Layer order matches
        // volume load order: index 0 = underlay, last = overlay (if present).
        let value = null;
        const vals = Array.isArray(data.values) ? data.values : [];
        if (vals.length > 0) {
            const pick = hasOverlay ? vals[vals.length - 1] : vals[0];
            if (pick && typeof pick.value === "number") value = pick.value;
        }

        let text = `MNI [${x} ${y} ${z}] mm   value: ${fmtValue(value)}`;

        // Region name from the atlas layer (its raw value is the integer code;
        // .value stays numeric even though the layer carries a label colormap).
        if (atlasInfo && atlasInfo.index >= 0 && vals.length > atlasInfo.index) {
            const av = vals[atlasInfo.index];
            const labels = atlasInfo.labels || [];
            if (av && typeof av.value === "number") {
                const code = Math.round(av.value);
                if (code >= 1 && code <= labels.length && labels[code - 1]) {
                    text += `   region: ${labels[code - 1]}`;
                }
                // Highlight only this region (no-op if the code is unchanged).
                if (atlasCtl) atlasCtl.onCode(code);
            }
        }

        out.textContent = text;
    };
}

/**
 * Boot a CANlab NiiVue viewer.
 *
 * @param {(string|HTMLCanvasElement)} target
 *        Canvas element id (default "gl") or a canvas DOM element.
 * @param {object} [config]
 *        Viewer configuration. All fields optional except an image source.
 * @param {(string|{base64:string,name?:string})} [config.underlay]
 *        Anatomical underlay: URL string or {base64,name} payload. If omitted,
 *        the overlay (if any) is shown alone.
 * @param {(string|{base64:string,name?:string})} [config.overlay]
 *        Statistical/mask overlay: URL string or {base64,name} payload.
 * @param {(string|{base64:string,name?:string})} [config.atlas]
 *        Integer label (atlas) volume: URL string or {base64,name} payload.
 *        Loaded as a hidden NIFTI_INTENT_LABEL layer that drives the crosshair
 *        region-name readout and the "Atlas regions" Off/Filled/Outline toggle.
 * @param {string[]} [config.atlasLabels]
 *        Region names for the atlas; atlasLabels[k] names integer code k+1.
 * @param {string} [config.atlasName]            Atlas name (informational).
 * @param {string} [config.colormap="inferno"]   Overlay positive colormap.
 * @param {string} [config.colormapNegative="winter"] Overlay negative colormap.
 * @param {number} [config.cal_min]              Overlay lower threshold.
 * @param {number} [config.cal_max]              Overlay upper threshold.
 * @param {number} [config.opacity=0.8]          Overlay opacity (semi-transparent by default).
 * @param {string} [config.sliceType="multiplanar"]
 *        multiplanar|axial|coronal|sagittal|render.
 * @param {boolean} [config.showColorbar=true]   Show the NiiVue colorbar.
 * @param {number} [config.crosshairWidth=0.5]   Crosshair line width.
 * @param {number} [config.textHeight=0.04]      Relative font size for colorbar/labels.
 * @param {number} [config.colorbarHeight=0.04]  Relative colorbar height.
 * @param {boolean} [config.showOrientationLabels=false] Show A/P/S/I/L/R letters.
 * @param {boolean} [config.showRenderInMultiplanar=true] Show the 3D render tile in multiplanar.
 * @param {number[]} [config.crosshairColor]     RGBA crosshair color.
 * @param {number[]} [config.backColor]          RGBA canvas background color.
 * @param {boolean} [config.controls=true]       Build the layout button bar + overlay tools.
 * @param {boolean} [config.allowLoadOverlay=true] Add a "Load overlay" file picker to the controls.
 * @param {boolean} [config.showOpacity=true]    Build the overlay/underlay opacity sliders.
 * @param {boolean} [config.showReadout=true]    Wire the crosshair readout.
 * @returns {Promise<Niivue>} the initialized NiiVue instance.
 */
export async function canlabNiivue(target, config) {
    const cfg = config || {};

    const sliceName = cfg.sliceType || "multiplanar";
    const showColorbar = cfg.showColorbar !== false;
    const useControls = cfg.controls !== false;
    const showReadout = cfg.showReadout !== false;
    const showOpacity = cfg.showOpacity !== false;

    // --- Construct the NiiVue instance -------------------------------------
    // isResizeCanvas defaults to true: NiiVue uses a ResizeObserver to keep the
    // WebGL drawing buffer matched to the canvas's CSS display size (and devicePixelRatio).
    // Do NOT disable it, or clicks map only to the unscaled top-left of the canvas.
    const nvOpts = {
        isColorbar: showColorbar,
        // Compact, report-friendly defaults:
        crosshairWidth: cfg.crosshairWidth !== undefined ? cfg.crosshairWidth : 0.5, // ~half the default (1)
        textHeight: cfg.textHeight !== undefined ? cfg.textHeight : 0.04,             // smaller colorbar/label font (default .06)
        colorbarHeight: cfg.colorbarHeight !== undefined ? cfg.colorbarHeight : 0.04, // slimmer colorbar (default .05)
        // Orientation letters (A/P/S/I/L/R) collide with the image at small sizes; off by default.
        isOrientationTextVisible: cfg.showOrientationLabels === true,
        // Always show the 3D render tile beside the slices in multiplanar view
        // (default AUTO hides it unless the frame is wide).
        multiplanarShowRender: cfg.showRenderInMultiplanar === false ? SHOW_RENDER.AUTO : SHOW_RENDER.ALWAYS
    };
    if (cfg.crosshairColor) nvOpts.crosshairColor = cfg.crosshairColor;
    if (cfg.backColor) nvOpts.backColor = cfg.backColor;

    const nv = new Niivue(nvOpts);

    // --- Attach to the canvas ----------------------------------------------
    if (typeof target === "string" || target === undefined) {
        await nv.attachTo(target || "gl");
    } else {
        nv.attachToCanvas(target);
    }

    // --- Load images: underlay first, then overlay -------------------------
    const hasUnderlay = !!cfg.underlay;
    const hasOverlay = !!cfg.overlay;
    // Default to a semi-transparent overlay so the anatomy reads through the blobs.
    const overlayOpacity = cfg.opacity !== undefined && cfg.opacity !== null ? cfg.opacity : 0.8;

    if (hasUnderlay) {
        // Underlay rendered in grayscale by default.
        await loadImage(nv, cfg.underlay, { colormap: "gray" }, "underlay.nii.gz");
    }

    // --- Atlas label volume (hidden; drives the region readout + show toggle) ---
    // Loaded BEFORE the overlay so the overlay stays the LAST volume — every
    // "last == overlay" assumption below (colormap dropdown, threshold slider,
    // load-overlay swap) keeps holding. The atlas starts hidden (opacity 0) but
    // NiiVue still samples its integer label value at the crosshair.
    const hasAtlas = !!cfg.atlas;
    let atlasIndex = -1;
    if (hasAtlas) {
        try {
            await loadImage(nv, cfg.atlas, { opacity: 0, colorbarVisible: false }, "atlas.nii.gz");
            atlasIndex = nv.volumes.length - 1;
            const av = nv.volumes[atlasIndex];
            if (av) {
                // Mark as a label map: NiiVue's integer-label shader, per-region
                // colormapLabel coloring, and the atlasOutline toggle all key off
                // NIFTI_INTENT_LABEL (intent_code 1002).
                if (av.hdr) av.hdr.intent_code = 1002;
                av.colorbarVisible = false;
                av.setColormapLabel(buildAtlasLabelColormap(cfg.atlasLabels || []));
                // Thicker outline band (read by the CANlab patch in niivue.js;
                // ignored -> ~1-voxel outline if that patch is ever lost).
                nv.opts.atlasOutlineWidth = ATLAS_OUTLINE_WIDTH;
                nv.setOpacity(atlasIndex, 0);
                nv.updateGLVolume();
            }
        } catch (e) {
            // A failed atlas layer must not break the viewer; the readout simply
            // omits the region name and the toggle button is not built.
            // eslint-disable-next-line no-console
            console.error("canlab_niivue: failed to set up atlas layer", e);
            atlasIndex = -1;
        }
    }

    let overlayIndex = -1;
    let overlayCalMin = NaN;
    let overlayCalMax = NaN;
    if (hasOverlay) {
        const overlayOpts = {
            colormap: cfg.colormap || "inferno",
            colormapNegative: cfg.colormapNegative || "winter",
            opacity: overlayOpacity
        };
        if (cfg.cal_min !== undefined && cfg.cal_min !== null) overlayOpts.cal_min = cfg.cal_min;
        if (cfg.cal_max !== undefined && cfg.cal_max !== null) overlayOpts.cal_max = cfg.cal_max;

        await loadImage(nv, cfg.overlay, overlayOpts, "overlay.nii.gz");

        // Make subthreshold (|value| < cal_min) and zero/NaN voxels transparent so
        // the overlay shows blobs only and the underlay shows through elsewhere.
        overlayIndex = nv.volumes.length - 1;
        const ov = nv.volumes[overlayIndex];
        if (ov) {
            ov.colormapType = COLORMAP_TYPE_TRANSPARENT_BELOW_MIN;
            if (overlayOpts.cal_min !== undefined) ov.cal_min = overlayOpts.cal_min;
            if (overlayOpts.cal_max !== undefined) ov.cal_max = overlayOpts.cal_max;
            nv.updateGLVolume();
            // Record the effective range for the threshold slider (caller-provided
            // values, or NiiVue's auto-computed cal_min/cal_max from the data).
            overlayCalMin = typeof ov.cal_min === "number" ? ov.cal_min : NaN;
            overlayCalMax = typeof ov.cal_max === "number" ? ov.cal_max : NaN;
        }
    }

    // --- Slice layout ------------------------------------------------------
    nv.setSliceType(sliceTypeFromName(sliceName));

    // --- Controls, opacity sliders, readout --------------------------------
    if (useControls) {
        buildControls(nv, sliceName);
        if (hasOverlay) buildOverlayControls(nv, cfg);
    }
    if (showOpacity) {
        buildOpacityPanel(nv, {
            hasUnderlay, hasOverlay, overlayIndex, overlayOpacity,
            calMin: overlayCalMin, calMax: overlayCalMax
        });
    }

    // Atlas-region control (dropdown + single-region highlighter), built after
    // the slider panel so it appears alongside the sliders.
    let atlasCtl = null;
    if (atlasIndex >= 0) {
        atlasCtl = buildAtlasControls(nv, atlasIndex);
    }

    if (showReadout) {
        wireReadout(nv, hasOverlay, { index: atlasIndex, labels: cfg.atlasLabels || [] }, atlasCtl);
    }

    // Seed the highlight from the initial crosshair position (NiiVue may not
    // fire onLocationChange until the first interaction).
    if (atlasCtl && atlasIndex >= 0) {
        try {
            const av = nv.volumes[atlasIndex];
            const mm = nv.frac2mm(nv.scene.crosshairPos);
            const vox = av.mm2vox(mm);
            const code = Math.round(av.getValue(vox[0], vox[1], vox[2]));
            atlasCtl.onCode(code);
        } catch (e) {
            // Non-fatal: the first crosshair move will highlight the region.
        }
    }

    return nv;
}
