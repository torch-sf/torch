"""
numerics/spectal
==========

Holds functions for generating a spectral field.
"""
__all__ = [
    "SpectralFieldGenerator"
]

import numpy as np

class SpectralFieldGenerator:
    """
    Generate random scalar or vector fields with a prescribed isotropic
    power-law spectrum.

    The class constructs a random spectral field (k-space), applies a
    band-limited power-law filter, and transforms the result to real space
    using an inverse FFT.

    The resulting field is statistically homogeneous and isotropic, with
    power spectrum:

        P(k) \propto k^slope

    where k = |k|.

    Note
    ----
    The generated fields are not divergence-free and do not enforce
    any physical constraints such as incompressibility.
    """

    def __init__(self, shape, kmin, kmax, slope, 
                 noise="uniform", seed = 0):
        """ 
        Initialized the class.
        
        Arguments
        ---------
            shape : tuple of int
                Grid size, e.g. (128, 128, 128).
            kmin : float
                Minimum wavenumber
            kmax : float
                Maximum wavenumber
            slope : float
                Power-law index of the power spectrum.
            noise : str
                Mode for noise in Fourier space (uniform or normal).
            seed : int
                Random number seed for sampling

        """
        self.shape = shape
        self.kmin = kmin
        self.kmax = kmax
        self.slope = slope
        self.noise = noise
        self.seed = seed
        self.rng = np.random.default_rng(seed)

    def _complex_noise(self):
        """
        Generate complex random noise in Fourier space.

        Supported modes:
            uniform: bounded white noise in [-0.5, 0.5]
            normal: normal distribution N(0,1)

        Returns
        -------
        field_k : ndarray (complex)
            Complex-valued Fourier-space field.
        """
        if self.noise == "uniform":
            return (self.rng.uniform(*self.shape) - 0.5) + 1j * (self.rng.uniform(*self.shape) - 0.5)
        elif self.noise == "normal":
            return (
                self.rng.normal(size=self.shape)
                + 1j * self.rng.normal(size=self.shape)
            ) / np.sqrt(2)

        else:
            raise ValueError(f"Unknown noise type: {self.noise}")

    def _compute_kgrid(self):
        """
        Compute the Fourier-space wavenumber grid.

        Constructs discrete wavevectors (kx, ky, kz) consistent with FFT
        and computes the wavenumbers:

            k^2 = kx^2 + ky^2 + kz^2

        Returns
        -------
        k : ndarray (float)
            Wavenumber field.
        """
        kx = np.fft.fftfreq(self.shape[0]) * self.shape[0]
        ky = np.fft.fftfreq(self.shape[1]) * self.shape[1]
        kz = np.fft.fftfreq(self.shape[2]) * self.shape[2]
        (kx, ky, kz) = np.meshgrid(kx,ky,kz, indexing="ij")
        return np.sqrt(kx**2 + ky**2 + kz**2)

    def _build_mask(self, k):
        """
        Construct band-limiting mask for wavenumbers kmin <= k <= kmax

        Arguments
        ----------
        k : ndarray
            Wavenumber field.

        Returns
        -------
        mask : ndarray (float)
            Real-valued spectral weighting function.
        """
        mask = np.ones_like(k)
        mask[(k < self.kmin) | (k > self.kmax)] = 0
        return mask

    def power_spectrum(self, k):
        """
        Apply amplitude scaling following power spectrum
            
            P(k) \propto k^slope,

        i.e., Fourier modes follow |\tilde{f}(k)| = \sqrt{P(k)}

        Arguments
        ----------
        k : ndarray
            Wavenumber field.

        Returns
        -------
        power : ndarray (float)
            Amplitudes rescaled to power spectrum.
        """
        k_safe = np.where(k == 0, 1, k)
        return k_safe**self.slope

    def generate_fourier_field(self):
        """
        Construct Fourier-space field with prescribed power spectrum.
        
        Returns
        -------
        field_k : ndarray (float)
            Fourier-space field.
        """
        k = self._compute_kgrid()

        mask = self._build_mask(k)
        amplitude = np.sqrt(self.power_spectrum(k))
        noise = self._complex_noise()
        
        field_k = noise * amplitude * mask
        return field_k
    
    def generate_scalar_field(self):
        """
        Generate a real-valued scalar field in configuration space.

        Returns
        -------
        field_s : ndarray (float)
            Real-valued scalar field
        """
        field_k = self.generate_fourier_field()
        return np.fft.ifftn(field_k).real
    
    def generate_vector_field(self, components=3):
        """
        Generate a multi-component vector field.
        Each component is independent and follows the same spectrum.

        Arguments
        ---------
        components : int
            Number of components

        Returns
        -------
        fields : tuple(ndarray)
            Multi-component vector field
        """
        fields = []

        for _ in range(components):
            field = self.generate_scalar_field()
            fields.append(field)

        return tuple(fields)
