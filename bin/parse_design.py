#!/usr/bin/env python3
import pandas as pd
import sys


df = pd.read_csv(sys.argv[1])
df.to_csv('design.csv', index=False)
