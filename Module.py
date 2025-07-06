import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    classification_report, confusion_matrix,
    roc_auc_score, RocCurveDisplay, roc_curve
)
from sklearn.utils import resample
from sklearn.preprocessing import label_binarize

from Bio import SeqIO
from collections import Counter
import math

import argparse
import gzip
import os
import shutil
from pathlib import Path
import marimo as mo

import gdown
import requests
