import { useState, useEffect, useCallback, useRef } from "react";
import type {
  ScatteringParams,
  ScatteringResult,
  FieldResult,
} from "../types/cylinder";

// Dynamic import for WASM module
let wasmModule: typeof import("scattering-core") | null = null;
let wasmInitPromise: Promise<void> | null = null;

async function initWasm(): Promise<typeof import("scattering-core")> {
  if (wasmModule) return wasmModule;

  if (!wasmInitPromise) {
    wasmInitPromise = (async () => {
      const module = await import("scattering-core");
      await module.default();
      wasmModule = module;
    })();
  }

  await wasmInitPromise;
  return wasmModule!;
}

export interface UseScatteringResult {
  isLoading: boolean;
  isReady: boolean;
  error: string | null;
  scatteringResult: ScatteringResult | null;
  fieldResult: FieldResult | null;
  computeScattering: (params: ScatteringParams) => Promise<ScatteringResult>;
  computeField: (
    params: ScatteringParams,
    scattering: ScatteringResult,
  ) => Promise<FieldResult>;
  computeAll: (
    params: ScatteringParams,
  ) => Promise<{ scattering: ScatteringResult; field: FieldResult }>;
}

export function useScattering(): UseScatteringResult {
  const [isLoading, setIsLoading] = useState(true);
  const [isReady, setIsReady] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [scatteringResult, setScatteringResult] =
    useState<ScatteringResult | null>(null);
  const [fieldResult, setFieldResult] = useState<FieldResult | null>(null);
  const moduleRef = useRef<typeof import("scattering-core") | null>(null);

  // Initialize WASM module
  useEffect(() => {
    let mounted = true;

    initWasm()
      .then((module) => {
        if (mounted) {
          moduleRef.current = module;
          setIsReady(true);
          setIsLoading(false);
          console.log("WASM module loaded:", module.get_info());
        }
      })
      .catch((err) => {
        if (mounted) {
          setError(`Failed to load WASM module: ${err.message}`);
          setIsLoading(false);
        }
      });

    return () => {
      mounted = false;
    };
  }, []);

  const computeScattering = useCallback(
    async (params: ScatteringParams): Promise<ScatteringResult> => {
      if (!moduleRef.current) {
        throw new Error("WASM module not loaded");
      }

      const t0 = performance.now();
      const result = moduleRef.current.compute_scattering(
        params.wavelength,
        params.material.permittivity.re,
        params.material.permittivity.im,
        params.material.permeability.re,
        params.material.permeability.im,
        params.polarization,
        params.maxOrder,
      );
      console.log(
        `Scattering computation: ${(performance.now() - t0).toFixed(1)}ms`,
      );

      setScatteringResult(result);
      return result;
    },
    [],
  );

  const computeField = useCallback(
    async (
      params: ScatteringParams,
      scattering: ScatteringResult,
    ): Promise<FieldResult> => {
      if (!moduleRef.current) {
        throw new Error("WASM module not loaded");
      }

      const numCoeffs = scattering.orders.length;

      // Convert coefficients to regular arrays (WASM expects Vec<f64> and Vec<i32>)
      const incidentReal: number[] = [];
      const incidentImag: number[] = [];
      const scatteringReal: number[] = [];
      const scatteringImag: number[] = [];
      const internalReal: number[] = [];
      const internalImag: number[] = [];
      const orders: number[] = [];

      for (let i = 0; i < numCoeffs; i++) {
        // Complex numbers are [re, im] tuples from WASM
        const inc = scattering.incident_coefficients[i];
        const scat = scattering.scattering_coefficients[i];
        const int = scattering.internal_coefficients[i];

        incidentReal.push(inc[0]);
        incidentImag.push(inc[1]);
        scatteringReal.push(scat[0]);
        scatteringImag.push(scat[1]);
        internalReal.push(int[0]);
        internalImag.push(int[1]);
        orders.push(scattering.orders[i]);
      }

      const t0 = performance.now();
      const result = moduleRef.current.compute_electric_field(
        params.wavelength,
        params.material.permittivity.re,
        params.material.permittivity.im,
        params.material.permeability.re,
        params.material.permeability.im,
        new Float64Array(incidentReal),
        new Float64Array(incidentImag),
        new Float64Array(scatteringReal),
        new Float64Array(scatteringImag),
        new Float64Array(internalReal),
        new Float64Array(internalImag),
        new Int32Array(orders),
      );
      console.log(
        `Field computation: ${(performance.now() - t0).toFixed(1)}ms`,
      );

      setFieldResult(result);
      return result;
    },
    [],
  );

  const computeAll = useCallback(
    async (params: ScatteringParams) => {
      const scattering = await computeScattering(params);
      const field = await computeField(params, scattering);
      return { scattering, field };
    },
    [computeScattering, computeField],
  );

  return {
    isLoading,
    isReady,
    error,
    scatteringResult,
    fieldResult,
    computeScattering,
    computeField,
    computeAll,
  };
}
