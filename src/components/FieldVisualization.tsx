import type { JSX } from "react";
import { useRef, useEffect, useState, useMemo } from "react";
import type { FieldResult, Polarization } from "../types/cylinder";
import "./FieldVisualization.css";

export type VisualizationMode = "magnitude" | "real" | "imag" | "phase";

interface FieldVisualizationProps {
  fieldResult: FieldResult | null;
  isComputing: boolean;
  polarization: Polarization;
  width?: number;
  height?: number;
}

// Colormap functions
function getViridisColor(t: number): [number, number, number] {
  // Simplified viridis colormap
  t = Math.max(0, Math.min(1, t));
  const r = Math.round(
    255 *
      (0.267004 +
        t * (0.329415 + t * (-0.814464 + t * (2.28653 - t * 1.06868)))),
  );
  const g = Math.round(
    255 *
      (0.004874 +
        t * (0.873449 + t * (0.107514 + t * (-0.631923 + t * 0.645732)))),
  );
  const b = Math.round(
    255 *
      (0.329415 +
        t * (1.01541 + t * (-1.67917 + t * (1.59456 - t * 0.594634)))),
  );
  return [
    Math.max(0, Math.min(255, r)),
    Math.max(0, Math.min(255, g)),
    Math.max(0, Math.min(255, b)),
  ];
}

function getPhaseColor(phase: number): [number, number, number] {
  // HSL to RGB for cyclic phase colormap
  const h = ((phase + Math.PI) / (2 * Math.PI)) * 360;
  const s = 0.8;
  const l = 0.5;

  const c = (1 - Math.abs(2 * l - 1)) * s;
  const x = c * (1 - Math.abs(((h / 60) % 2) - 1));
  const m = l - c / 2;

  let r = 0,
    g = 0,
    b = 0;
  if (h < 60) {
    r = c;
    g = x;
    b = 0;
  } else if (h < 120) {
    r = x;
    g = c;
    b = 0;
  } else if (h < 180) {
    r = 0;
    g = c;
    b = x;
  } else if (h < 240) {
    r = 0;
    g = x;
    b = c;
  } else if (h < 300) {
    r = x;
    g = 0;
    b = c;
  } else {
    r = c;
    g = 0;
    b = x;
  }

  return [
    Math.round((r + m) * 255),
    Math.round((g + m) * 255),
    Math.round((b + m) * 255),
  ];
}

function getDivergingColor(t: number): [number, number, number] {
  // Blue-white-red diverging colormap for real/imag parts
  t = Math.max(0, Math.min(1, t));
  if (t < 0.5) {
    const s = t * 2;
    return [
      Math.round(59 + s * (255 - 59)),
      Math.round(76 + s * (255 - 76)),
      Math.round(192 + s * (255 - 192)),
    ];
  } else {
    const s = (t - 0.5) * 2;
    return [
      Math.round(255 - s * (255 - 180)),
      Math.round(255 - s * (255 - 59)),
      Math.round(255 - s * (255 - 59)),
    ];
  }
}

export function FieldVisualization({
  fieldResult,
  isComputing,
  polarization,
  width = 512,
  height = 512,
}: FieldVisualizationProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [mode, setMode] = useState<VisualizationMode>("magnitude");

  // Compute visualization values and stats using useMemo
  const { values, stats } = useMemo(() => {
    if (!fieldResult) {
      return { values: null, stats: null };
    }

    const { field_real, field_imag, grid_size } = fieldResult;
    const vals = new Float64Array(grid_size * grid_size);
    let minVal = Infinity;
    let maxVal = -Infinity;

    for (let i = 0; i < vals.length; i++) {
      const re = field_real[i];
      const im = field_imag[i];

      let val: number;
      switch (mode) {
        case "magnitude":
          val = Math.sqrt(re * re + im * im);
          break;
        case "real":
          val = re;
          break;
        case "imag":
          val = im;
          break;
        case "phase":
          val = Math.atan2(im, re);
          break;
      }
      vals[i] = val;

      if (mode !== "phase") {
        minVal = Math.min(minVal, val);
        maxVal = Math.max(maxVal, val);
      }
    }

    if (mode === "phase") {
      minVal = -Math.PI;
      maxVal = Math.PI;
    }

    return { values: vals, stats: { min: minVal, max: maxVal } };
  }, [fieldResult, mode]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !fieldResult || !values || !stats) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const { grid_size } = fieldResult;
    const { min: minVal, max: maxVal } = stats;

    // Create image data
    const imageData = ctx.createImageData(grid_size, grid_size);
    const data = imageData.data;

    for (let i = 0; i < values.length; i++) {
      let t: number;
      let color: [number, number, number];

      if (mode === "phase") {
        color = getPhaseColor(values[i]);
      } else if (mode === "real" || mode === "imag") {
        // Diverging colormap centered on zero
        const absMax = Math.max(Math.abs(minVal), Math.abs(maxVal));
        t = absMax > 0 ? (values[i] + absMax) / (2 * absMax) : 0.5;
        color = getDivergingColor(t);
      } else {
        // Magnitude: sequential colormap
        t = maxVal > minVal ? (values[i] - minVal) / (maxVal - minVal) : 0;
        color = getViridisColor(t);
      }

      const idx = i * 4;
      data[idx] = color[0];
      data[idx + 1] = color[1];
      data[idx + 2] = color[2];
      data[idx + 3] = 255;
    }

    // Draw to a temporary canvas at full resolution
    const tempCanvas = document.createElement("canvas");
    tempCanvas.width = grid_size;
    tempCanvas.height = grid_size;
    const tempCtx = tempCanvas.getContext("2d")!;
    tempCtx.putImageData(imageData, 0, 0);

    // Scale to display canvas
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(tempCanvas, 0, 0, width, height);

    // Draw cylinder outline
    const centerX = width / 2;
    const centerY = height / 2;
    const cylinderRadius = (0.5 / fieldResult.view_size) * width;

    ctx.strokeStyle = "rgba(255, 255, 255, 0.8)";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(centerX, centerY, cylinderRadius, 0, 2 * Math.PI);
    ctx.stroke();

    // Draw axis markers
    ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
    ctx.lineWidth = 1;
    ctx.setLineDash([5, 5]);

    // Horizontal axis
    ctx.beginPath();
    ctx.moveTo(0, centerY);
    ctx.lineTo(width, centerY);
    ctx.stroke();

    // Vertical axis
    ctx.beginPath();
    ctx.moveTo(centerX, 0);
    ctx.lineTo(centerX, height);
    ctx.stroke();

    ctx.setLineDash([]);
  }, [fieldResult, mode, width, height, values, stats]);

  // Use E_z for TM polarization, H_z for TE polarization
  const fieldSymbol = polarization === "TM" ? "E" : "H";
  const modeLabels: Record<VisualizationMode, JSX.Element> = {
    magnitude: (
      <span>
        |{fieldSymbol}
        <sub>z</sub>|
      </span>
    ),
    real: (
      <span>
        Re({fieldSymbol}
        <sub>z</sub>)
      </span>
    ),
    imag: (
      <span>
        Im({fieldSymbol}
        <sub>z</sub>)
      </span>
    ),
    phase: (
      <span>
        Phase({fieldSymbol}
        <sub>z</sub>)
      </span>
    ),
  };

  return (
    <div className="field-visualization">
      <div className="visualization-controls">
        <label>Display:</label>
        <div className="mode-buttons">
          {(Object.keys(modeLabels) as VisualizationMode[]).map((m) => (
            <button
              key={m}
              className={mode === m ? "active" : ""}
              onClick={() => setMode(m)}
            >
              {modeLabels[m]}
            </button>
          ))}
        </div>
      </div>

      <div className="canvas-container">
        {isComputing && (
          <div className="computing-overlay">
            <div className="spinner" />
            <span>Computing field...</span>
          </div>
        )}
        <canvas
          ref={canvasRef}
          width={width}
          height={height}
          style={{ opacity: isComputing ? 0.5 : 1 }}
        />
      </div>

      {stats && fieldResult && (
        <div className="field-stats">
          <span>
            Grid: {fieldResult.grid_size}×{fieldResult.grid_size}
          </span>
          <span>
            Range: [{stats.min.toExponential(2)}, {stats.max.toExponential(2)}]
          </span>
        </div>
      )}

      {!fieldResult && !isComputing && (
        <div className="no-data">
          <p>Click "Compute Field" to calculate the electromagnetic field</p>
        </div>
      )}
    </div>
  );
}
