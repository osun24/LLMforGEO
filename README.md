# Installation
## Create virtual environment (recommended)
python -m venv llm_geo_env
source llm_geo_env/bin/activate  # On Windows: llm_geo_env\Scripts\activate

## Install dependencies
pip install -r requirements.txt

# Requirements
- Python 3.8+ (developed using 3.10.9)

## Dependencies
- pandas>=2.0.0
- beautifulsoup4>=4.12.0
- biopython>=1.80
- openai>=1.50.0
- selenium>=4.20.0 (relies on ChromeDriver to retrieve clinical data from GEO2R)
- tiktoken>=0.8.0
- openpyxl>=3.1.0
- python-dotenv>=1.0.0
- PySide6>=6.8.0
- requests>=2.32.0
