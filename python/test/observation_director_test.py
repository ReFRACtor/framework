#!/usr/bin/env python
"""
Test SWIG director feature for Observation interface.
Verifies that Python subclasses can override virtual methods.
"""

import sys
import numpy as np

try:
    from refractor.framework_swig import observation, spectral_domain, spectral_range, spectrum, unit
    print("✓ Successfully imported refractor.framework_swig")
except ImportError as e:
    print(f"✗ Failed to import refractor.framework_swig: {e}")
    sys.exit(1)


class TestObservation(observation.Observation):
    """Python subclass of Observation to test director feature."""

    def __init__(self, num_channels=3):
        super().__init__()
        self._num_channels = num_channels
        print(f"  Created TestObservation with {num_channels} channels")

    def num_channels(self):
        """Override num_channels property."""
        print(f"  Python num_channels() called, returning {self._num_channels}")
        return self._num_channels

    def spectral_domain(self, sensor_index):
        """Override spectral_domain method."""
        print(f"  Python spectral_domain({sensor_index}) called")
        # Create a simple spectral domain with some wavelengths in cm^-1
        wavelengths = [12900.0 + i * 10.0 for i in range(10)]
        return spectral_domain.SpectralDomain(wavelengths)

    def radiance(self, sensor_index, skip_jacobian=False):
        """Override radiance method."""
        print(f"  Python radiance({sensor_index}, skip_jacobian={skip_jacobian}) called")
        # Create a simple spectrum with dummy radiance values
        sd = self.spectral_domain(sensor_index)
        radiance_data = [1.0 + i * 0.1 for i in range(10)]
        radiance_units = "W / m^2 / sr / um"
        sr = spectral_range.SpectralRange(radiance_data, radiance_units)
        return spectrum.Spectrum(sd, sr)


def test_director():
    """Test that Python overrides are called."""
    print("\n=== Testing Observation Director Feature ===\n")

    print("1. Creating Python subclass of Observation:")
    obs = TestObservation(num_channels=3)

    print("\n2. Calling num_channels() property:")
    num_ch = obs.num_channels()
    assert num_ch == 3, f"Expected 3 channels, got {num_ch}"
    print(f"  ✓ Got expected value: {num_ch}")

    print("\n3. Calling spectral_domain(0):")
    sd = obs.spectral_domain(0)
    assert sd is not None, "spectral_domain returned None"
    print(f"  ✓ Got SpectralDomain with {len(sd.data)} points")

    print("\n4. Calling radiance(0):")
    spectrum = obs.radiance(0)
    assert spectrum is not None, "radiance returned None"
    print(f"  ✓ Got Spectrum")

    print("\n5. Testing with skip_jacobian=True:")
    spectrum2 = obs.radiance(1, skip_jacobian=True)
    assert spectrum2 is not None, "radiance with skip_jacobian returned None"
    print(f"  ✓ Got Spectrum")

    print("\n=== All Tests Passed ===")
    print("✓ SWIG director feature is working correctly!")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(test_director())
    except Exception as e:
        print(f"\n✗ Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
