"""Test module."""
from pathlib import Path

tests_path = Path(__file__).resolve().parents[0]
golden_data = Path(tests_path, 'golden_data')