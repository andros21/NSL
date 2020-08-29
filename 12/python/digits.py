#!/usr/bin/env python
# -*- coding: utf-8 -*-

# MANIPULATION OF GIMP HAND-WRITTEN DIGITS

import os
import numpy as np
from PIL import Image

# List of GIMP manipulated images
gimp_digits = []

for digit in os.listdir('digits/'):
    
    digit_in = Image.open(                   # Load gimp handwritten digits,
            f'digits/{digit}').convert('1')  # and convert to bitmap
    
    digit_np = np.asarray(digit_in)          # From Image to Array
    
    digit_np = np.invert(digit_np)           # Invert boolean
    
    gimp_digits.append(digit_np)
