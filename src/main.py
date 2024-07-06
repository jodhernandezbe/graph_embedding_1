#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Main module.

This module contains the main function to run the complete pipeline.

"""

from src.data_processing import run_processing


def run() -> None:
    """Run the complete pipeline."""
    run_processing()


if __name__ == "__main__":
    run()
