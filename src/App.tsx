import { useState, useEffect } from "react";
import { CylinderControls } from "./components/CylinderControls";
import { FieldVisualization } from "./components/FieldVisualization";
import { useScattering } from "./hooks/useScattering";
import {
  createDefaultParams,
  calculateSizeParameter,
  calculateRefractiveIndex,
  formatComplex,
} from "./types/cylinder";
import type { ScatteringParams } from "./types/cylinder";
import "./App.css";

function App() {
  const [params, setParams] = useState<ScatteringParams>(createDefaultParams());
  const [isComputing, setIsComputing] = useState(false);

  const {
    isLoading,
    isReady,
    error,
    scatteringResult,
    fieldResult,
    computeAll,
  } = useScattering();

  // Auto-compute when params change or when WASM becomes ready
  useEffect(() => {
    if (!isReady) return;

    // Computation is fast enough (~2-5ms) to run immediately
    const compute = async () => {
      setIsComputing(true);
      try {
        await computeAll(params);
      } catch (err) {
        console.error("Computation error:", err);
      } finally {
        setIsComputing(false);
      }
    };

    compute();
  }, [isReady, params, computeAll]);

  return (
    <div className="app">
      <header className="app-header">
        <h1>2D Electromagnetic Scattering</h1>
        <p>
          Interactive visualization of EM wave scattering from an infinite
          cylinder
        </p>
      </header>

      <main className="app-main">
        <div className="visualization-panel">
          <FieldVisualization
            fieldResult={fieldResult}
            isComputing={isComputing}
            polarization={params.polarization}
            width={512}
            height={512}
          />

          <div className="status-section">
            {isLoading && (
              <span className="status">Loading WASM module...</span>
            )}
            {isComputing && <span className="status">Computing...</span>}
            {error && <span className="status error">{error}</span>}
          </div>

          <div className="scattering-info">
            <h3>Scattering Info</h3>
            <p>
              Size parameter: x ={" "}
              {calculateSizeParameter(params.wavelength).toFixed(4)}
            </p>
            <p>
              Permittivity: εᵣ = {formatComplex(params.material.permittivity)}
            </p>
            <p>
              Permeability: μᵣ = {formatComplex(params.material.permeability)}
            </p>
            <p>
              Refractive index: n ={" "}
              {formatComplex(calculateRefractiveIndex(params.material))}
            </p>
            {scatteringResult && (
              <p>
                Orders: {scatteringResult.orders[0]} to{" "}
                {scatteringResult.orders[scatteringResult.orders.length - 1]} (
                {scatteringResult.orders.length} terms)
              </p>
            )}
          </div>
        </div>

        <div className="controls-panel">
          <CylinderControls params={params} onChange={setParams} />
        </div>
      </main>
    </div>
  );
}

export default App;
