import { useRef, useEffect } from "react";
import type { ScatteringParams } from "../types/cylinder";

interface ScatteringCanvasProps {
  params: ScatteringParams;
  width?: number;
  height?: number;
}

// Fixed cylinder diameter
const CYLINDER_DIAMETER = 1.0;
const CYLINDER_RADIUS = CYLINDER_DIAMETER / 2;

/**
 * Canvas component that renders the 2D xy-plane with the cylinder cross-section.
 */
export function ScatteringCanvas({
  params,
  width = 500,
  height = 500,
}: ScatteringCanvasProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    // Clear canvas
    ctx.fillStyle = "#1a1a2e";
    ctx.fillRect(0, 0, width, height);

    // Calculate scale: show view extending beyond the cylinder
    const viewRadius = 1.5; // Show ±1.5 units (cylinder has radius 0.5)
    const scale = Math.min(width, height) / (2 * viewRadius);

    // Transform to center origin
    const centerX = width / 2;
    const centerY = height / 2;

    // Draw grid
    drawGrid(ctx, centerX, centerY, scale, viewRadius);

    // Draw incident wave direction indicator
    drawIncidentWave(ctx, centerX, centerY, scale, viewRadius, params);

    // Draw axes
    drawAxes(ctx, centerX, centerY, width, height);

    // Draw cylinder
    drawCylinder(ctx, centerX, centerY, scale, params);

    // Draw labels
    drawLabels(ctx, centerX, centerY, scale, viewRadius);
  }, [params, width, height]);

  return (
    <canvas
      ref={canvasRef}
      width={width}
      height={height}
      style={{ border: "1px solid #333", borderRadius: "4px" }}
    />
  );
}

function drawGrid(
  ctx: CanvasRenderingContext2D,
  centerX: number,
  centerY: number,
  scale: number,
  viewRadius: number,
) {
  ctx.strokeStyle = "#2a2a4a";
  ctx.lineWidth = 1;

  const gridSpacing = 0.5;

  // Draw vertical lines
  for (let x = -viewRadius; x <= viewRadius; x += gridSpacing) {
    const screenX = centerX + x * scale;
    ctx.beginPath();
    ctx.moveTo(screenX, 0);
    ctx.lineTo(screenX, ctx.canvas.height);
    ctx.stroke();
  }

  // Draw horizontal lines
  for (let y = -viewRadius; y <= viewRadius; y += gridSpacing) {
    const screenY = centerY - y * scale;
    ctx.beginPath();
    ctx.moveTo(0, screenY);
    ctx.lineTo(ctx.canvas.width, screenY);
    ctx.stroke();
  }
}

function drawAxes(
  ctx: CanvasRenderingContext2D,
  centerX: number,
  centerY: number,
  width: number,
  height: number,
) {
  ctx.strokeStyle = "#666";
  ctx.lineWidth = 2;

  // X-axis
  ctx.beginPath();
  ctx.moveTo(0, centerY);
  ctx.lineTo(width, centerY);
  ctx.stroke();

  // Y-axis
  ctx.beginPath();
  ctx.moveTo(centerX, 0);
  ctx.lineTo(centerX, height);
  ctx.stroke();

  // Axis labels
  ctx.fillStyle = "#888";
  ctx.font = "14px monospace";
  ctx.fillText("x", width - 20, centerY - 10);
  ctx.fillText("y", centerX + 10, 20);
}

function drawIncidentWave(
  ctx: CanvasRenderingContext2D,
  centerX: number,
  centerY: number,
  scale: number,
  viewRadius: number,
  params: ScatteringParams,
) {
  // Draw arrow indicating incident wave direction (coming from -x direction)
  const arrowX = -viewRadius + 0.2;
  const arrowLength = 0.4;

  ctx.strokeStyle = "#ffaa00";
  ctx.fillStyle = "#ffaa00";
  ctx.lineWidth = 2;

  const startX = centerX + arrowX * scale;
  const endX = centerX + (arrowX + arrowLength) * scale;
  const y = centerY;

  // Arrow line
  ctx.beginPath();
  ctx.moveTo(startX, y);
  ctx.lineTo(endX, y);
  ctx.stroke();

  // Arrow head
  ctx.beginPath();
  ctx.moveTo(endX, y);
  ctx.lineTo(endX - 10, y - 6);
  ctx.lineTo(endX - 10, y + 6);
  ctx.closePath();
  ctx.fill();

  // Label
  ctx.font = "11px monospace";
  ctx.fillText("k", startX + 5, y - 10);

  // Polarization indicator
  const polLabel = params.polarization === "TM" ? "E∥z" : "H∥z";
  ctx.fillText(polLabel, startX, y + 20);
}

function drawCylinder(
  ctx: CanvasRenderingContext2D,
  centerX: number,
  centerY: number,
  scale: number,
  params: ScatteringParams,
) {
  const screenRadius = CYLINDER_RADIUS * scale;

  // Fill cylinder with color based on material properties
  const { permittivity } = params.material;

  // Color intensity based on real part of permittivity
  const alpha = Math.min(0.6, 0.1 + Math.abs(permittivity.re) * 0.05);

  // Hue: blue for positive permittivity, orange for negative (like metals)
  const hue = permittivity.re >= 0 ? 220 : 30;

  ctx.fillStyle = `hsla(${hue}, 70%, 50%, ${alpha})`;
  ctx.beginPath();
  ctx.arc(centerX, centerY, screenRadius, 0, 2 * Math.PI);
  ctx.fill();

  // Draw cylinder outline
  ctx.strokeStyle = "#00d4ff";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(centerX, centerY, screenRadius, 0, 2 * Math.PI);
  ctx.stroke();

  // Draw center point
  ctx.fillStyle = "#ff6b6b";
  ctx.beginPath();
  ctx.arc(centerX, centerY, 4, 0, 2 * Math.PI);
  ctx.fill();

  // Draw diameter indicator
  ctx.strokeStyle = "#ff6b6b";
  ctx.lineWidth = 1;
  ctx.setLineDash([5, 5]);
  ctx.beginPath();
  ctx.moveTo(centerX - screenRadius, centerY);
  ctx.lineTo(centerX + screenRadius, centerY);
  ctx.stroke();
  ctx.setLineDash([]);

  // Diameter label
  ctx.fillStyle = "#ff6b6b";
  ctx.font = "12px monospace";
  ctx.fillText("d = 1", centerX - 15, centerY + screenRadius + 20);
}

function drawLabels(
  ctx: CanvasRenderingContext2D,
  centerX: number,
  centerY: number,
  scale: number,
  viewRadius: number,
) {
  const gridSpacing = 0.5;

  ctx.fillStyle = "#666";
  ctx.font = "10px monospace";

  // X-axis labels
  for (let x = -viewRadius; x <= viewRadius; x += gridSpacing) {
    if (Math.abs(x) < 0.001) continue;
    const screenX = centerX + x * scale;
    ctx.fillText(x.toFixed(1), screenX - 10, centerY + 15);
  }

  // Y-axis labels
  for (let y = -viewRadius; y <= viewRadius; y += gridSpacing) {
    if (Math.abs(y) < 0.001) continue;
    const screenY = centerY - y * scale;
    ctx.fillText(y.toFixed(1), centerX + 5, screenY + 4);
  }
}
