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
 * @returns {void}
 */
function wireReadout(nv, hasOverlay) {
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

        out.textContent = `MNI [${x} ${y} ${z}] mm   value: ${fmtValue(value)}`;
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
    if (showReadout) wireReadout(nv, hasOverlay);

    return nv;
}
