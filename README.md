# TCSAC: Enhanced Point Cloud Registration for Workpieces Using Triangular Constraints in Complex Industrial Scenarios 🌐🧠

[![Python 3.6+](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-Paper-<COLOR>.svg)](https://arxiv.org/abs/XXXX.XXXX)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yourusername/yourrepo/)

> **Three-Dimensional Correspondence Selection for Automated Calibration**  
> Advanced point cloud registration for 3D matching and industrial workpiece alignment

## 📌 Table of Contents
- [Features](#-features)
- [Installation](#-installation)
- [Usage](#-usage)
- [Dataset](#-dataset)
- [Configuration](#-configuration)
- [Calibration](#-calibration)
- [Testing](#-testing)
- [Citation](#-citation)
- [License](#-license)
- [Contact](#-contact)

## ✨ Features
- 🎯 **High-Precision Registration**: Robust algorithm for 3D point cloud alignment
- ⚡ **Optimized Performance**: Efficient correspondence selection
- 🛠️ **Industrial Ready**: Specialized workpiece configuration
- 📊 **Comprehensive Logging**: Detailed test execution records
- 🤖 **Robot-Camera Calibration**: Sphere-based calibration method

## 🛠️ Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/TCSAC.git
   cd TCSAC
2. Install dependencies:
   ```bash
    pip install -r requirements.txt
## 🚀 Usage
1. Main Algorithm Execution:
   ```bash
    python TCSAC.py \
        --config config_json/3DMatch_config.json \
        --input test/3DMatch_sample \
        --output results/3DMatch_output
## 📂 Dataset
1. The `test/` directory contains partial test data. Complete datasets are available upon request:
    ```text
        test/
        ├── 3DMatch/      # Sample 3DMatch point clouds
        │   ├── fragments
        │   └── gt_result
        └── workpiece/    # Sample industrial workpieces
            ├── fragments
            └── gt_result
## ⚙️ Configuration
1. Customize parameters in JSON config files:
   ```bash
    // config_json/3DMatch_config.json
    {
        "feature_radius": 0.3,
        "max_correspondence": 500,
        "ransac_iterations": 10000,
        "confidence_threshold": 0.99
    }
## 🔍 Calibration
1. For robot-camera calibration:
   ```bash
    cd sphere-calibration
    python calibrate.py --images calibration_images/ --output calibration.json
## 🧪 Testing
1. Run test scripts with logging:
   ```bash
    # For 3DMatch evaluation
    python TCSAC_test_3DMatch.py > logs/3DMatch_test.log 2>&1
    # For workpiece evaluation
    python TCSAC_test_workpiece.py > logs/workpiece_test.log 2>&1
## 📄 License
This project is licensed under the MIT License - see the LICENSE file for details.

