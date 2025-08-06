import os, re, time, json, sys, subprocess
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup
import openai, tiktoken
from dotenv import load_dotenv
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                               QPushButton, QTextEdit, QLineEdit, QLabel, QFormLayout, QGroupBox, QMessageBox, QProgressBar, QSpacerItem, QSizePolicy, QCheckBox, QDialog, QDialogButtonBox, QSplitter)  # NEW: Added QSplitter
from PySide6.QtCore import QThread, Signal, QSettings, Qt, QMetaObject, Q_ARG, Qt, Slot  # NEW: Added Qt
from PySide6.QtGui import QRegularExpressionValidator
from PySide6.QtCore import QRegularExpression

# Global settings
settings = QSettings("MyCompany", "NCBIApp")
Entrez.email = settings.value("email", "", type=str)
openai.api_key = settings.value("openai_api_key", "", type=str)

# Global OpenAI client (set later in handle_batch)
client = None

# Global callback for logging to the UI (set in MainWindow)
global_log_callback = None

FIXED_STUDIES = ["samples_table_GSE29013.csv","samples_table_GSE6044.csv","samples_table_GSE14814.csv","samples_table_GSE7880.csv","samples_table_GSE39279.csv","samples_table_GSE50081.csv","samples_table_GSE39345.csv","samples_table_GSE42127.csv","samples_table_GSE37745.csv","samples_table_GSE47115.csv","samples_table_GSE77209.csv","samples_table_GSE61926.csv","samples_table_GSE30291.csv","samples_table_GSE13255.csv","samples_table_GSE84146.csv","samples_table_GSE189047.csv","samples_table_GSE115246.csv","samples_table_GSE10096.csv","samples_table_GSE27606.csv","samples_table_GSE6914.csv","samples_table_GSE38054.csv","samples_table_GSE27389.csv","samples_table_GSE74706.csv","samples_table_GSE240758.csv","samples_table_GSE8332.csv","samples_table_GSE107246.csv","samples_table_GSE57083.csv","samples_table_GSE9212.csv","samples_table_GSE229301.csv","samples_table_GSE84797.csv","samples_table_GSE92679.csv","samples_table_GSE75960.csv","samples_table_GSE36681.csv","samples_table_GSE11969.csv","samples_table_GSE55859.csv","samples_table_GSE115864.csv","samples_table_GSE121323.csv","samples_table_GSE27388.csv","samples_table_GSE52717.csv","samples_table_GSE18385.csv","samples_table_GSE171702.csv","samples_table_GSE21276.csv","samples_table_GSE18346.csv","samples_table_GSE68851.csv","samples_table_GSE42076.csv","samples_table_GSE157062.csv","samples_table_GSE108214.csv","samples_table_GSE26241.csv","samples_table_GSE62359.csv","samples_table_GSE145880.csv","samples_table_GSE16194.csv","samples_table_GSE32036.csv","samples_table_GSE40517.csv","samples_table_GSE136932.csv","samples_table_GSE184414.csv","samples_table_GSE114761.csv","samples_table_GSE31852.csv","samples_table_GSE229302.csv","samples_table_GSE115457.csv","samples_table_GSE203514.csv","samples_table_GSE75037.csv","samples_table_GSE115456.csv","samples_table_GSE161116.csv","samples_table_GSE23066.csv","samples_table_GSE26644.csv","samples_table_GSE136933.csv","samples_table_GSE118370.csv","samples_table_GSE183178.csv","samples_table_GSE171517.csv","samples_table_GSE32989.csv","samples_table_GSE216508.csv","samples_table_GSE50659.csv","samples_table_GSE249568.csv","samples_table_GSE49644.csv","samples_table_GSE189045.csv","samples_table_GSE28282.csv","samples_table_GSE222124.csv","samples_table_GSE24876.csv","samples_table_GSE10245.csv","samples_table_GSE9994.csv","samples_table_GSE197109.csv","samples_table_GSE175601.csv","samples_table_GSE37765.csv","samples_table_GSE69631.csv","samples_table_GSE37759.csv","samples_table_GSE40738.csv","samples_table_GSE142278.csv","samples_table_GSE54981.csv","samples_table_GSE43580.csv","samples_table_GSE209983.csv","samples_table_GSE203510.csv","samples_table_GSE33072.csv","samples_table_GSE77803.csv","samples_table_GSE121090.csv","samples_table_GSE47056.csv","samples_table_GSE64472.csv","samples_table_GSE69630.csv","samples_table_GSE40275.csv","samples_table_GSE69181.csv","samples_table_GSE127191.csv","samples_table_GSE87879.csv","samples_table_GSE21933.csv","samples_table_GSE106765.csv","samples_table_GSE54162.csv","samples_table_GSE59.csv","samples_table_GSE121682.csv","samples_table_GSE17073.csv","samples_table_GSE68869.csv","samples_table_GSE168280.csv","samples_table_GSE188406.csv","samples_table_GSE33845.csv","samples_table_GSE94393.csv","samples_table_GSE165029.csv","samples_table_GSE189268.csv","samples_table_GSE50138.csv","samples_table_GSE43622.csv","samples_table_GSE44617.csv","samples_table_GSE136934.csv","samples_table_GSE16025.csv","samples_table_GSE248249.csv","samples_table_GSE103780.csv","samples_table_GSE11945.csv","samples_table_GSE43568.csv","samples_table_GSE202156.csv","samples_table_GSE131518.csv","samples_table_GSE6695.csv","samples_table_GSE179212.csv","samples_table_GSE16597.csv","samples_table_GSE63685.csv","samples_table_GSE72788.csv","samples_table_GSE176444.csv","samples_table_GSE55860.csv","samples_table_GSE198238.csv","samples_table_GSE136935.csv","samples_table_GSE69747.csv","samples_table_GSE52712.csv","samples_table_GSE6253.csv","samples_table_GSE58571.csv","samples_table_GSE31935.csv","samples_table_GSE31946.csv","samples_table_GSE95766.csv","samples_table_GSE22047.csv","samples_table_GSE154286.csv","samples_table_GSE7670.csv","samples_table_GSE166750.csv","samples_table_GSE75309.csv","samples_table_GSE62113.csv","samples_table_GSE98929.csv","samples_table_GSE117181.csv","samples_table_GSE101085.csv","samples_table_GSE44077.csv","samples_table_GSE37700.csv","samples_table_GSE103527.csv","samples_table_GSE27705.csv","samples_table_GSE37271.csv","samples_table_GSE58661.csv","samples_table_GSE131027.csv","samples_table_GSE79308.csv","samples_table_GSE137481.csv","samples_table_GSE42375.csv","samples_table_GSE75308.csv","samples_table_GSE103888.csv","samples_table_GSE57156.csv","samples_table_GSE11117.csv","samples_table_GSE46729.csv","samples_table_GSE10547.csv","samples_table_GSE64007.csv","samples_table_GSE66604.csv","samples_table_GSE86928.csv","samples_table_GSE193719.csv","samples_table_GSE70919.csv","samples_table_GSE28582.csv","samples_table_GSE93300.csv","samples_table_GSE66606.csv","samples_table_GSE123031.csv","samples_table_GSE119144.csv","samples_table_GSE14079.csv","samples_table_GSE158695.csv","samples_table_GSE33363.csv","samples_table_GSE30118.csv","samples_table_GSE31428.csv","samples_table_GSE27262.csv","samples_table_GSE84720.csv","samples_table_GSE69482.csv","samples_table_GSE23206.csv","samples_table_GSE116699.csv","samples_table_GSE78217.csv","samples_table_GSE21612.csv","samples_table_GSE37058.csv","samples_table_GSE63384.csv","samples_table_GSE79486.csv","samples_table_GSE67727.csv","samples_table_GSE37138.csv","samples_table_GSE8569.csv","samples_table_GSE43493.csv","samples_table_GSE4186.csv","samples_table_GSE51722.csv","samples_table_GSE1987.csv","samples_table_GSE31798.csv","samples_table_GSE154243.csv","samples_table_GSE64766.csv","samples_table_GSE42979.csv","samples_table_GSE37699.csv","samples_table_GSE147029.csv","samples_table_GSE242843.csv","samples_table_GSE160482.csv","samples_table_GSE142623.csv","samples_table_GSE196499.csv","samples_table_GSE67051.csv","samples_table_GSE101929.csv","samples_table_GSE69732.csv","samples_table_GSE3202.csv","samples_table_GSE137445.csv","samples_table_GSE137479.csv","samples_table_GSE22874.csv","samples_table_GSE54712.csv","samples_table_GSE32497.csv","samples_table_GSE112376.csv","samples_table_GSE94536.csv","samples_table_GSE122538.csv","samples_table_GSE84094.csv","samples_table_GSE207715.csv","samples_table_GSE84095.csv","samples_table_GSE62504.csv","samples_table_GSE269782.csv","samples_table_GSE161400.csv","samples_table_GSE16534.csv","samples_table_GSE32496.csv","samples_table_GSE71587.csv","samples_table_GSE60887.csv","samples_table_GSE101684.csv","samples_table_GSE108139.csv","samples_table_GSE42373.csv","samples_table_GSE75468.csv","samples_table_GSE110825.csv","samples_table_GSE7878.csv","samples_table_GSE4342.csv","samples_table_GSE134788.csv","samples_table_GSE188665.csv","samples_table_GSE66616.csv","samples_table_GSE31799.csv","samples_table_GSE80344.csv","samples_table_GSE42749.csv","samples_table_GSE57422.csv","samples_table_GSE33198.csv","samples_table_GSE12630.csv","samples_table_GSE43494.csv","samples_table_GSE118274.csv","samples_table_GSE49422.csv","samples_table_GSE27489.csv","samples_table_GSE166999.csv","samples_table_GSE11729.csv","samples_table_GSE76759.csv","samples_table_GSE29249.csv","samples_table_GSE101862.csv","samples_table_GSE164750.csv","samples_table_GSE5843.csv","samples_table_GSE99993.csv","samples_table_GSE22863.csv","samples_table_GSE24584.csv","samples_table_GSE112375.csv","samples_table_GSE84040.csv","samples_table_GSE78210.csv","samples_table_GSE19188.csv","samples_table_GSE29077.csv","samples_table_GSE72715.csv","samples_table_GSE84096.csv","samples_table_GSE112374.csv","samples_table_GSE22862.csv","samples_table_GSE15002.csv","samples_table_GSE57966.csv","samples_table_GSE101863.csv","samples_table_GSE29248.csv","samples_table_GSE123340.csv","samples_table_GSE129692.csv","samples_table_GSE166998.csv","samples_table_GSE20318.csv","samples_table_GSE138172.csv","samples_table_GSE40407.csv","samples_table_GSE23361.csv","samples_table_GSE9315.csv","samples_table_GSE112214.csv","samples_table_GSE83666.csv","samples_table_GSE197236.csv","samples_table_GSE193707.csv","samples_table_GSE71358.csv","samples_table_GSE228854.csv","samples_table_GSE75466.csv","samples_table_GSE25251.csv","samples_table_GSE131952.csv","samples_table_GSE34228.csv","samples_table_GSE43249.csv","samples_table_GSE228049.csv","samples_table_GSE103512.csv","samples_table_GSE116679.csv","samples_table_GSE127472.csv","samples_table_GSE27902.csv","samples_table_GSE73162.csv","samples_table_GSE54293.csv","samples_table_GSE29250.csv","samples_table_GSE108136.csv","samples_table_GSE147153.csv","samples_table_GSE65933.csv","samples_table_GSE28827.csv","samples_table_GSE125113.csv","samples_table_GSE75467.csv","samples_table_GSE122126.csv","samples_table_GSE54495.csv","samples_table_GSE71403.csv","samples_table_GSE37318.csv","samples_table_GSE1037.csv","samples_table_GSE41271.csv","samples_table_GSE77925.csv","samples_table_GSE43458.csv","samples_table_GSE53893.csv","samples_table_GSE71398.csv","samples_table_GSE213816.csv","samples_table_GSE114694.csv","samples_table_GSE42791.csv","samples_table_GSE4824.csv","samples_table_GSE74948.csv","samples_table_GSE25508.csv","samples_table_GSE14925.csv","samples_table_GSE19434.csv","samples_table_GSE6410.csv","samples_table_GSE136780.csv","samples_table_GSE23019.csv","samples_table_GSE25326.csv","samples_table_GSE73160.csv","samples_table_GSE137475.csv","samples_table_GSE134381.csv","samples_table_GSE45175.csv","samples_table_GSE25118.csv","samples_table_GSE19804.csv","samples_table_GSE108492.csv","samples_table_GSE26438.csv","samples_table_GSE51854.csv","samples_table_GSE81644.csv","samples_table_GSE73161.csv","samples_table_GSE38404.csv","samples_table_GSE57000.csv","samples_table_GSE108135.csv","samples_table_GSE42425.csv","samples_table_GSE52144.csv","samples_table_GSE166997.csv","samples_table_GSE110815.csv","samples_table_GSE196494.csv","samples_table_GSE28993.csv","samples_table_GSE198446.csv","samples_table_GSE12428.csv","samples_table_GSE50627.csv","samples_table_GSE119400.csv","samples_table_GSE99143.csv","samples_table_GSE43459.csv","samples_table_GSE53882.csv","samples_table_GSE56044.csv","samples_table_GSE28571.csv","samples_table_GSE29135.csv","samples_table_GSE51501.csv","samples_table_GSE35911.csv","samples_table_GSE141479.csv","samples_table_GSE65854.csv","samples_table_GSE40424.csv","samples_table_GSE83490.csv","samples_table_GSE19592.csv","samples_table_GSE59997.csv","samples_table_GSE101673.csv","samples_table_GSE181067.csv","samples_table_GSE201608.csv","samples_table_GSE102178.csv","samples_table_GSE100753.csv","samples_table_GSE27911.csv","samples_table_GSE73414.csv","samples_table_GSE31800.csv","samples_table_GSE36821.csv","samples_table_GSE81641.csv","samples_table_GSE73158.csv","samples_table_GSE3593.csv","samples_table_GSE204942.csv","samples_table_GSE159431.csv","samples_table_GSE87026.csv","samples_table_GSE27284.csv","samples_table_GSE36216.csv","samples_table_GSE31579.csv","samples_table_GSE75307.csv","samples_table_GSE51266.csv","samples_table_GSE19034.csv","samples_table_GSE5519.csv","samples_table_GSE38944.csv","samples_table_GSE13525.csv","samples_table_GSE20853.csv","samples_table_GSE155157.csv","samples_table_GSE42998.csv","samples_table_GSE103891.csv","samples_table_GSE35912.csv","samples_table_GSE135304.csv","samples_table_GSE20304.csv","samples_table_GSE67729.csv","samples_table_GSE45403.csv","samples_table_GSE5123.csv","samples_table_GSE52143.csv","samples_table_GSE225959.csv","samples_table_GSE81643.csv","samples_table_GSE68793.csv","samples_table_GSE104244.csv","samples_table_GSE127462.csv","samples_table_GSE31625.csv","samples_table_GSE51853.csv","samples_table_GSE51852.csv","samples_table_GSE189607.csv","samples_table_GSE60480.csv","samples_table_GSE81642.csv","samples_table_GSE35640.csv","samples_table_GSE73167.csv","samples_table_GSE180347.csv","samples_table_GSE18842.csv","samples_table_GSE67728.csv","samples_table_GSE2922.csv","samples_table_GSE31552.csv","samples_table_GSE2088.csv","samples_table_GSE166749.csv","samples_table_GSE249262.csv","samples_table_GSE48244.csv","samples_table_GSE72197.csv","samples_table_GSE69574.csv","samples_table_GSE189270.csv","samples_table_GSE29391.csv","samples_table_GSE93157.csv","samples_table_GSE136131.csv","samples_table_GSE101979.csv","samples_table_GSE263003.csv","samples_table_GSE106151.csv","samples_table_GSE81143.csv","samples_table_GSE16559.csv","samples_table_GSE92843.csv","samples_table_GSE15723.csv","samples_table_GSE47059.csv","samples_table_GSE4716.csv","samples_table_GSE9971.csv","samples_table_GSE79404.csv","samples_table_GSE7035.csv","samples_table_GSE30979.csv","samples_table_GSE165540.csv","samples_table_GSE70540.csv","samples_table_GSE69561.csv","samples_table_GSE25019.csv","samples_table_GSE56036.csv","samples_table_GSE80316.csv","samples_table_GSE72194.csv","samples_table_GSE69563.csv","samples_table_GSE147708.csv","samples_table_GSE272045.csv","samples_table_GSE26704.csv","samples_table_GSE79228.csv","samples_table_GSE223395.csv","samples_table_GSE13191.csv","samples_table_GSE54033.csv","samples_table_GSE43567.csv","samples_table_GSE126043.csv","samples_table_GSE244643.csv","samples_table_GSE93586.csv","samples_table_GSE24502.csv","samples_table_GSE48433.csv","samples_table_GSE33532.csv","samples_table_GSE163913.csv","samples_table_GSE138650.csv","samples_table_GSE94809.csv","samples_table_GSE58542.csv","samples_table_GSE10089.csv","samples_table_GSE123066.csv","samples_table_GSE130740.csv","samples_table_GSE15126.csv","samples_table_GSE146601.csv","samples_table_GSE72195.csv","samples_table_GSE66872.csv","samples_table_GSE77764.csv","samples_table_GSE133715.csv","samples_table_GSE98979.csv","samples_table_GSE102286.csv","samples_table_GSE50913.csv","samples_table_GSE66245.csv","samples_table_GSE5816.csv","samples_table_GSE198958.csv","samples_table_GSE244647.csv","samples_table_GSE244646.csv","samples_table_GSE84799.csv","samples_table_GSE168466.csv","samples_table_GSE198959.csv","samples_table_GSE8894.csv","samples_table_GSE52308.csv","samples_table_GSE252795.csv","samples_table_GSE37779.csv","samples_table_GSE66244.csv","samples_table_GSE101836.csv","samples_table_GSE38310.csv","samples_table_GSE16760.csv","samples_table_GSE102287.csv","samples_table_GSE264186.csv","samples_table_GSE12472.csv","samples_table_GSE2003.csv","samples_table_GSE51946.csv","samples_table_GSE99870.csv","samples_table_GSE4882.csv","samples_table_GSE143018.csv","samples_table_GSE63074.csv","samples_table_GSE4869.csv","samples_table_GSE104151.csv","samples_table_GSE138092.csv","samples_table_GSE65968.csv","samples_table_GSE169587.csv","samples_table_GSE104757.csv","samples_table_GSE84200.csv","samples_table_GSE27554.csv","samples_table_GSE33247.csv","samples_table_GSE66087.csv","samples_table_GSE115458.csv","samples_table_GSE17681.csv","samples_table_GSE126045.csv","samples_table_GSE63882.csv","samples_table_GSE244645.csv","samples_table_GSE193628.csv","samples_table_GSE16575.csv","samples_table_GSE68929.csv","samples_table_GSE158940.csv","samples_table_GSE63459.csv","samples_table_GSE64322.csv","samples_table_GSE21656.csv","samples_table_GSE55869.csv","samples_table_GSE5828.csv","samples_table_GSE117702.csv","samples_table_GSE190731.csv","samples_table_GSE51945.csv","samples_table_GSE116959.csv","samples_table_GSE44390.csv","samples_table_GSE61676.csv","samples_table_GSE5579.csv"]

def log_message(message, verbose_only=False):
    is_verbose = settings.value("verbose", False, type=bool)
    if verbose_only and not is_verbose:
        return
    print(message)
    if global_log_callback:
        # Use the bound object's metaObject to ensure the slot is called in the main thread.
        QMetaObject.invokeMethod(global_log_callback.__self__, 'append_log', Qt.QueuedConnection, Q_ARG(str, message))


class SettingsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        layout = QFormLayout(self)

        self.email_edit = QLineEdit()
        self.api_key_edit = QLineEdit()
        self.api_key_edit.setEchoMode(QLineEdit.Password)
        # NEW: Highlight API key textbox in red if not set
        current_api = settings.value("openai_api_key", "", type=str)
        if not current_api.strip():
            self.api_key_edit.setStyleSheet("border: 2px solid red;")
        else:
            self.api_key_edit.setStyleSheet("")
        # NEW: Add user-input for max entries to fetch
        self.retmax_edit = QLineEdit()
        # NEW: Add verbose mode toggle
        self.verbose_checkbox = QCheckBox()
        
        # Load saved settings
        self.email_edit.setText(settings.value("email", "", type=str))
        self.api_key_edit.setText(settings.value("openai_api_key", "", type=str))
        self.retmax_edit.setText(str(settings.value("retmax", "10000", type=str)))  # default 10000
        self.verbose_checkbox.setChecked(settings.value("verbose", False, type=bool))
        
        layout.addRow("Email:", self.email_edit)
        layout.addRow("API Key:", self.api_key_edit)
        layout.addRow("Max Entries to Fetch:", self.retmax_edit)
        layout.addRow("Verbose Logging:", self.verbose_checkbox)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.save_settings)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def save_settings(self):
        settings.setValue("email", self.email_edit.text())
        settings.setValue("openai_api_key", self.api_key_edit.text())
        settings.setValue("retmax", self.retmax_edit.text())  # NEW
        settings.setValue("verbose", self.verbose_checkbox.isChecked())
        # Update global values
        Entrez.email = self.email_edit.text()
        openai.api_key = self.api_key_edit.text()
        self.accept()

def save_metadata_to_excel(metadata_list, filename="geo_metadata.xlsx"):
    metadata_df = pd.DataFrame(metadata_list)
    try:
        existing_df = pd.read_excel(filename, engine='openpyxl')
        metadata_df = pd.concat([existing_df, metadata_df], ignore_index=True)
    except FileNotFoundError:
        pass
    metadata_df.to_excel(filename, index=False, engine='openpyxl')
    log_message(f"Metadata saved to {filename}")

def search_geo(term, retmax=5):
    handle = Entrez.esearch(db="gds", term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_geo_metadata(geo_id, retries=3, delay=0.3):
    for attempt in range(retries):
        try:
            time.sleep(delay)
            handle = Entrez.efetch(db="gds", id=geo_id, rettype="xml")
            metadata = handle.read()
            handle.close()
            return metadata
        except Exception as e:
            log_message(f"Error fetching metadata for GEO ID {geo_id} (Attempt {attempt+1}/{retries}): {e}", verbose_only=True)
            time.sleep(2 ** attempt)
    log_message(f"Failed to fetch metadata for GEO ID {geo_id} after {retries} retries.")
    return None

def extract_description(metadata):
    if pd.isnull(metadata): return None
    split_point = metadata.find("Organism:")
    return metadata[:split_point].strip() if split_point != -1 else metadata.strip()

def extract_field(metadata, field_name):
    if pd.isnull(metadata): return None
    field_marker = f"{field_name}:"
    if field_marker in metadata:
        try:
            field_start = metadata.index(field_marker) + len(field_marker)
            remaining = metadata[field_start:]
            next_field = re.search(r"\b\w+:", remaining)
            return metadata[field_start: next_field.start() + field_start].strip() if next_field else metadata[field_start:].strip()
        except ValueError:
            return None
    return None

def extract_samples(metadata):
    if pd.isnull(metadata): return None
    match = re.search(r'(\d+)\s+Samples', metadata)
    return int(match.group(1)) if match else None

def clean_platforms(metadata):
    if pd.isnull(metadata): return None
    return re.sub(r'\s+\d+\s+Samples', '', metadata).strip()

def process_final_excel(input_file, output_file):
    data = pd.read_excel(input_file)
    fields = {
        'Organism': 'Organism',
        'Type': 'Type',
        'Platform(s)': ['Platform', 'Platforms'],
        'FTP download': 'FTP download',
        'SeriesAccession': ['Series Accession', 'Accession'],
        'ID': 'ID'
    }
    data['Description'] = data['Metadata'].apply(extract_description)
    for col, field in fields.items():
        if isinstance(field, list):
            data[col] = data['Metadata'].apply(lambda x: extract_field(x, field[0]) or extract_field(x, field[1]))
        else:
            data[col] = data['Metadata'].apply(lambda x: extract_field(x, field))
    data['Samples'] = data['Platform(s)'].apply(extract_samples)
    data['Platform(s)'] = data['Platform(s)'].apply(clean_platforms)
    data.to_excel(output_file, index=False)
    log_message(f"Processed data saved to {output_file}")

def run_ncbi_search(search_query, retmax, excel_file, final_file, progress_callback=None):
    geo_ids = search_geo(search_query, retmax=retmax)
    log_message(f"Number of GEO datasets found: {len(geo_ids)}")
    if not geo_ids:
        log_message("No GEO datasets found.")
        return
    metadata_list = []
    total_ids = len(geo_ids)
    for idx, geo_id in enumerate(geo_ids):
        log_message(f"Fetching metadata for GEO ID {geo_id} ({idx+1}/{total_ids})", verbose_only=True)
        metadata = fetch_geo_metadata(geo_id)
        if metadata:
            log_message(f"Metadata fetched for GEO ID {geo_id}")
            metadata_list.append({"GEO ID": geo_id, "Metadata": metadata})
        if (idx + 1) % 10 == 0 or idx == total_ids - 1:
            log_message(f"Saving metadata for {idx+1} GEO datasets...", verbose_only=True)
            save_metadata_to_excel(metadata_list, filename=excel_file)
            metadata_list = []
        if (idx + 1) % 100 == 0:
            log_message("Pausing to prevent rate-limiting...", verbose_only=True)
            time.sleep(5)
        if progress_callback:
            progress = int((idx+1) / total_ids * 100)
            progress_callback(progress, f"Fetched {idx+1}/{total_ids} datasets")
    process_final_excel(excel_file, final_file)
    if progress_callback:
        progress_callback(100, "NCBI data fetch completed.")

# === excel-to-geo2r-csv Functions ===
error_log_file = "error_log.txt"
no_data = []

def collect_data_for_url(driver, accession_id, platform):
    if accession_id in no_data:
        log_message(f"Ignoring {accession_id} as it has no data.", verbose_only=True)
        return None
    url = (f"https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc={accession_id}&platform={platform}"
           if platform else f"https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc={accession_id}")
    log_message(f"Navigating to URL: {url}", verbose_only=True)
    driver.get(url)
    driver.implicitly_wait(10)
    try:
        WebDriverWait(driver, 10).until(EC.alert_is_present())
        alert = driver.switch_to.alert
        alert_text = alert.text
        alert.accept()
        with open(error_log_file, 'a') as log_file:
            log_file.write(f"Error for {accession_id} with platform {platform}: {alert_text}\n")
        return None
    except TimeoutException:
        pass
    html_content = driver.page_source
    soup = BeautifulSoup(html_content, 'html.parser')
    table_container = soup.find('div', class_='samplesTableContainer')
    if table_container:
        table = table_container.find('table')
        if table:
            rows = table.find_all('tr')
            headers, data_rows = [], []
            for i, row in enumerate(rows):
                if i == 0:
                    headers = [th.get_text(strip=True) for th in row.find_all('th')]
                else:
                    cells = [td.get_text(strip=True) for td in row.find_all('td')]
                    if cells:
                        data_rows.append(cells)
            return pd.DataFrame(data_rows, columns=headers)
    return None

def process_geo_data(input_file, output_folder, progress_callback=None):
    output_folder = output_folder or os.path.splitext(input_file)[0] + "_outputs"
    if not os.path.exists(input_file):
        log_message(f"Input file {input_file} not found. Creating an empty file with required columns.")
        empty_df = pd.DataFrame(columns=['SeriesAccession', 'Platform(s)'])
        empty_df.to_excel(input_file, index=False, engine='openpyxl')
    data = pd.read_excel(input_file, sheet_name='Sheet1')
    os.makedirs(output_folder, exist_ok=True)
    no_platform_file = os.path.join(output_folder, 'no_platforms.txt')
    driver = webdriver.Chrome()  # Ensure ChromeDriver is installed and in PATH
    no_data_file = os.path.join(output_folder, 'no_data.txt')
    global no_data
    no_data = []
    if os.path.exists(no_platform_file):
        with open(no_platform_file, 'r') as np_file:
            no_data = [line.split(':')[0] for line in np_file.readlines()]
    if os.path.exists(no_data_file):
        with open(no_data_file, 'r') as nd_file:
            no_data += [line.split(':')[0] for line in nd_file.readlines()]
    total_rows = len(data)
    for count, (_, row) in enumerate(data.iterrows(), start=1):
        accession_id = row['SeriesAccession']
        platforms = str(row['Platform(s)']).strip()
        output_file = os.path.join(output_folder, f'samples_table_{accession_id}.csv')
        log_message(f"Processing {accession_id} ({count}/{total_rows})", verbose_only=True)
        if os.path.exists(output_file) or accession_id in no_data:
            log_message(f"Skipping {accession_id} as it was processed already.", verbose_only=True)
            continue
        if platforms.lower() == 'nan' or not platforms:
            with open(no_platform_file, 'a') as np_file:
                np_file.write(f"{accession_id}: No platforms available.\n")
            continue
        platform_list = platforms.split()
        platform_dfs = []
        for platform in platform_list:
            df = collect_data_for_url(driver, accession_id, platform)
            if df is not None:
                df.columns = [f"{platform}_{col}" for col in df.columns]
                platform_dfs.append(df)
        if platform_dfs:
            combined_df = pd.concat(platform_dfs, axis=1)
            combined_df.to_csv(output_file, index=False, encoding='utf-8')
            log_message(f"Data successfully saved to {output_file}")
        else:
            log_message(f"No data found for {accession_id} with platforms: {platforms}")
            with open(no_data_file, 'a') as nd_file:
                nd_file.write(f"{accession_id}: No data found for platforms: {platforms}\n")
            time.sleep(1)
        if progress_callback:
            progress = int(count / total_rows * 100)
            progress_callback(progress, f"Processed {count}/{total_rows} rows")
    driver.quit()
    if progress_callback:
        progress_callback(100, "GEO data processing completed.")

def run_process_geo_data(input_file, output_folder, progress_callback=None):
    process_geo_data(input_file, output_folder, progress_callback)

# === csv-to-4o-mini Functions ===
def count_tokens(messages, model="gpt-4o-mini"):
    encoding = tiktoken.encoding_for_model(model)
    total_tokens = 0
    for message in messages:
        total_tokens += len(encoding.encode(message["role"])) + len(encoding.encode(message["content"]))
    return total_tokens

def create_batch_file(input_directory, batch_file_path, system_message, user_message):
    tasks = []
    total_tokens = 0
    for filename in os.listdir(input_directory):
        if filename.endswith(".csv"):
            csv_path = os.path.join(input_directory, filename)
            df = pd.read_csv(csv_path)
            csv_text = json.dumps(df.to_csv(index=False))
            combined_user_message = f"{user_message}\n\n{csv_text}"
            messages = [
                {"role": "system", "content": system_message},
                {"role": "user", "content": combined_user_message}
            ]
            task = {
                "custom_id": filename,
                "method": "POST",
                "url": "/v1/chat/completions",
                "body": {"model": "gpt-4o-mini", "messages": messages, "temperature": 0}
            }
            task_tokens = count_tokens(messages, model="gpt-4o-mini")
            if task_tokens > 128000:
                log_message(f"{filename}: {task_tokens} tokens. Skipping because it exceeds the token limit.")
                continue
            tasks.append(task)
            total_tokens += task_tokens
            log_message(f"{filename}: {task_tokens} tokens", verbose_only=True)
    log_message(f"\nTotal estimated tokens for all requests: {total_tokens}")
    with open(batch_file_path, 'w') as batch_file:
        for task in tasks:
            batch_file.write(json.dumps(task) + '\n')
    return total_tokens

def upload_batch_file(batch_file_path):
    file = client.files.create(file=open(batch_file_path, 'rb'), purpose="batch")
    return file.id

def create_batch_job(file_id):
    response = client.batches.create(
        input_file_id=file_id,
        endpoint="/v1/chat/completions",
        completion_window="24h"
    )
    return response.id

def download_file(file_id, save_path):
    response = client.files.retrieve(file_id)
    file = client.files.content(file_id).content
    with open(save_path, 'wb') as f:
        f.write(file)

def run_create_batch_job(system_message, user_message, input_directory, batch_file_path, output_prefix, progress_callback=None):
    total_tokens = create_batch_file(input_directory, batch_file_path, system_message, user_message)
    log_message(f"Total tokens required: {total_tokens}. Proceeding with batch job creation...")
    file_id = upload_batch_file(batch_file_path)
    log_message(f"Batch file uploaded. File ID: {file_id}")
    batch_id = create_batch_job(file_id)
    log_message(f"Batch job ID: {batch_id}")
    log_message("Monitoring batch job status...")

    # Save last batch job details in settings
    settings.setValue("last_batch_id", batch_id)
    settings.setValue("last_output_prefix", output_prefix)
    
    max_wait = 1200  # 20 minutes
    elapsed = 0
    while True:
        response = client.batches.retrieve(batch_id)
        status = response.status
        # CHANGE THIS LATER
        if progress_callback:
            progress = min(100, int((elapsed / max_wait) * 100))
            progress_callback(progress, f"Waiting... {elapsed}s elapsed")
        log_message(f"Current status: {status}", verbose_only=True)
        log_message(f"Errors: {response.errors}", verbose_only=True)
        if status == "completed":
            if response.output_file_id:
                download_file(response.output_file_id, f"{output_prefix}.jsonl")
                log_message("Batch processing completed. Results saved.")
            if response.error_file_id:
                download_file(response.error_file_id, f"{output_prefix}_errors.jsonl")
                log_message("Some errors occurred. Details saved.")
            break
        elif status in {"failed", "cancelled", "expired"}:
            if response.error_file_id:
                download_file(response.error_file_id, f"{output_prefix}_errors.jsonl")
                log_message("Batch job ended with errors. Details saved.")
            break
        time.sleep(60)
        elapsed += 60

# === 4o-mini-to-excel Functions ===
def parse_batch_output(file_path, csv_output):
    rows = []
    for line in open(file_path, 'r'):
        data = json.loads(line)
        content = data['response']['body']['choices'][0]['message']['content']
        model = data['response']['body']['model']
        answers = {}
        for num, ans in re.findall(r'(\d+)\.\s*(.+)', content):
            answers[f"Q{num}"] = ans.strip()
        row = {"custom_id": data.get('custom_id', ''), "error": data.get('error', ''), "model": model}
        row.update(answers)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(csv_output, index=False)
    log_message(f"Parsed batch output saved to {csv_output}")

def run_parse_batch_output(jsonl_file, csv_output):
    parse_batch_output(jsonl_file, csv_output)

# === Worker Thread to Run Tasks Without Blocking the UI ===
class WorkerThread(QThread):
    finished_signal = Signal(str)
    def __init__(self, func, *args, **kwargs):
        super().__init__()
        self.func = func
        self.args = args
        self.kwargs = kwargs  # May include progress_callback
    def run(self):
        try:
            self.func(*self.args, **self.kwargs)
            self.finished_signal.emit("Task completed successfully.")
        except Exception as e:
            self.finished_signal.emit(f"Task failed: {e}")

# === Main Window Using PySide6 ===
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LLM Study Screening: NCBI GEO")
        central = QWidget()
        self.setCentralWidget(central)
        
        # --- Use QSplitter for adjustable heights ---
        main_layout = QVBoxLayout(central)
        splitter = QSplitter(Qt.Vertical)  # NEW: Allows user to adjust heights
        
        # Top container: project settings, NCBI and Batch settings, buttons
        top_widget = QWidget()
        top_layout = QVBoxLayout(top_widget)
        
        # --- Project Settings ---
        project_group = QGroupBox("Project Settings")
        project_layout = QHBoxLayout()
        project_name_label = QLabel("Project Name:")
        self.project_name_edit = QLineEdit("defaultProject")
        self.project_name_edit.textChanged.connect(self.handle_project_name_changed)
        project_layout.addWidget(project_name_label)
        project_layout.addWidget(self.project_name_edit)
        project_layout.addStretch()
        
        self.dark_mode_button = QPushButton("🌙")
        self.dark_mode_button.setFlat(True)
        self.dark_mode_button.setStyleSheet("""
            QPushButton {
                background: transparent;
                border: none;
                font-size: 16px;
                color: gray;
                padding: 4px;
            }
            QPushButton:hover { color: black; }
        """)
        self.dark_mode_button.clicked.connect(self.toggle_dark_mode)
        
        self.settings_button = QPushButton("⚙️")
        self.settings_button.setFlat(True)
        self.settings_button.setStyleSheet("""
            QPushButton {
                background: transparent;
                border: none;
                font-size: 16px;
                color: gray;
                padding: 4px;
            }
            QPushButton:hover { color: black; }
        """)
        self.settings_button.clicked.connect(self.open_settings_dialog)
        
        project_layout.addWidget(self.dark_mode_button)
        project_layout.addWidget(self.settings_button)
        project_group.setLayout(project_layout)
        top_layout.addWidget(project_group)
        
        # --- Menu Bar ---
        menu = self.menuBar().addMenu("File")
        settings_action = menu.addAction("Settings")
        settings_action.triggered.connect(self.open_settings_dialog)

        # --- NCBI Search Settings ---
        ncbi_group = QGroupBox("NCBI Search Settings")
        ncbi_group.setCheckable(True)
        ncbi_group.setChecked(True)

        ncbi_layout = QFormLayout()

        ncbi_group.toggled.connect(lambda checked, lay=ncbi_layout: 
            [lay.itemAt(i).widget().setVisible(checked) for i in range(lay.count())]
        )

        self.ncbi_query_edit = QTextEdit()
        # NEW: Allow user-adjustable height for text boxes
        self.ncbi_query_edit.setMinimumHeight(80)
        self.ncbi_query_edit.setMaximumHeight(100)
        self.ncbi_query_edit.setPlainText(
            '(("lung adenocarcinomas" OR "lung carcinomas" OR "lung squamous cell carcinomas" OR '
            '"large cell carcinomas" OR "bronchioalveolar carcinomas" OR "pulmonary adenocarcinomas" OR '
            '"lung adenocarcinoma" OR "lung carcinoma" OR "lung squamous cell carcinoma" OR "large cell carcinoma" OR '
            '"bronchioalveolar carcinoma" OR "pulmonary adenocarcinoma" OR "NSCLC" OR "non-small cell lung cancer") '
            'AND Homo sapiens[Organism] AND "gse"[Filter] AND ("expression profiling by array"[DataSet Type] OR '
            '"expression profiling by genome tiling array"[DataSet Type] OR "genome binding/occupancy profiling by array"[DataSet Type] OR '
            '"genome binding/occupancy profiling by genome tiling array"[DataSet Type] OR "genome variation profiling by array"[DataSet Type] OR '
            '"genome variation profiling by genome tiling array"[DataSet Type] OR "methylation profiling by array"[DataSet Type] OR '
            '"methylation profiling by genome tiling array"[DataSet Type] OR "non coding RNA profiling by array"[DataSet Type] OR '
            '"non coding RNA profiling by genome tiling array"[DataSet Type]))'
        )
        # NEW: Load retmax from settings for main window
        default_retmax = settings.value("retmax", "10000", type=str)
        self.ncbi_retmax_edit = QLineEdit(default_retmax)
        ncbi_layout.addRow(QLabel("Search Query:"), self.ncbi_query_edit)
        ncbi_layout.addRow(QLabel("Max Entries to Fetch:"), self.ncbi_retmax_edit)  # NEW
        ncbi_group.setLayout(ncbi_layout)
        top_layout.addWidget(ncbi_group)

        # --- Batch Job Settings ---
        batch_group = QGroupBox("Batch Job Settings")
        batch_layout = QFormLayout()
        self.batch_system_prompt_edit = QTextEdit()
        # NEW: Allow user-adjustable height for text boxes
        self.batch_system_prompt_edit.setMinimumHeight(50)
        self.batch_system_prompt_edit.setPlainText(
            "You are an oncology expert evaluating clinical datasets based on inclusion/exclusion criteria. Respond concisely and clearly, returning answers in the same numbered list format as the questions."
        )
        self.batch_user_prompt_edit = QTextEdit()
        # NEW: Allow user-adjustable height for text boxes
        self.batch_user_prompt_edit.setMinimumHeight(50)
        self.batch_user_prompt_edit.setPlainText(
            "Provide any additional context or description regarding the dataset."
        )
        batch_layout.addRow(QLabel("System Prompt:"), self.batch_system_prompt_edit)
        batch_layout.addRow(QLabel("User Prompt:"), self.batch_user_prompt_edit)
        
        # NEW: Container for individual model questions
        self.batch_questions_layout = QVBoxLayout()
        self.batch_questions_container = QWidget()
        self.batch_questions_container.setLayout(self.batch_questions_layout)
        batch_layout.addRow(QLabel("Queries:"), self.batch_questions_container)
        self.batch_question_edits = []
        self.add_question_button = QPushButton("Add Question")
        self.add_question_button.clicked.connect(self.add_question_field)
        batch_layout.addRow(self.add_question_button)
        self.add_question_field()
        
        batch_group.setLayout(batch_layout)
        top_layout.addWidget(batch_group)

        # --- Action Buttons ---
        button_layout = QHBoxLayout()
        self.ncbi_button = QPushButton("Fetch NCBI Data")
        self.geo_button = QPushButton("Process GEO Data")
        self.batch_button = QPushButton("Create Batch Job")
        self.parse_button = QPushButton("Get Last Batch")
        button_layout.addWidget(self.ncbi_button)
        button_layout.addWidget(self.geo_button)
        button_layout.addWidget(self.batch_button)
        button_layout.addWidget(self.parse_button)
        top_layout.addLayout(button_layout)
        
        # Add top widget to splitter
        splitter.addWidget(top_widget)
        
        # Bottom container: Log Output and Progress Bar
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout(bottom_widget)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        bottom_layout.addWidget(self.log_text)
        
        # --- New: Progress Bar ---
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("Progress: %p%")
        bottom_layout.addWidget(self.progress_bar)
        
        splitter.addWidget(bottom_widget)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        
        main_layout.addWidget(splitter)
        
        # Connect buttons to handlers
        self.ncbi_button.clicked.connect(self.handle_ncbi)
        self.geo_button.clicked.connect(self.handle_geo)
        self.batch_button.clicked.connect(self.handle_batch)
        self.parse_button.clicked.connect(self.handle_parse)
        
        # Load saved user input values from settings
        self.load_user_inputs()
        
        # Set the global log callback so that log_message() can append to the UI.
        global global_log_callback
        global_log_callback = self.append_log
        
        # Initialize dark_mode attribute
        self.dark_mode = is_macos_dark_mode() if sys.platform == "darwin" else False
        apply_styles(self.dark_mode)
    
    def load_user_inputs(self):
        s = QSettings("MyCompany", "NCBIApp")
        self.project_name_edit.setText(s.value("project_name", self.project_name_edit.text(), type=str))
        self.ncbi_query_edit.setPlainText(s.value("ncbi_query", self.ncbi_query_edit.toPlainText(), type=str))
        self.ncbi_retmax_edit.setText(s.value("ncbi_retmax", self.ncbi_retmax_edit.text(), type=str))
        self.batch_system_prompt_edit.setPlainText(s.value("batch_system_prompt", self.batch_system_prompt_edit.toPlainText(), type=str))
        self.batch_user_prompt_edit.setPlainText(s.value("batch_user_prompt", self.batch_user_prompt_edit.toPlainText(), type=str))
    
    def closeEvent(self, event):
        s = QSettings("MyCompany", "NCBIApp")
        s.setValue("project_name", self.project_name_edit.text())
        s.setValue("ncbi_query", self.ncbi_query_edit.toPlainText())
        s.setValue("ncbi_retmax", self.ncbi_retmax_edit.text())
        s.setValue("batch_system_prompt", self.batch_system_prompt_edit.toPlainText())
        s.setValue("batch_user_prompt", self.batch_user_prompt_edit.toPlainText())
        event.accept()
    
    def append_log(self, msg):
        self.log_text.append(msg)
    
    def update_progress(self, value, message):
        self.progress_bar.setValue(value)
        self.append_log(message)

    def handle_ncbi(self):
        project_prefix = self.project_name_edit.text()
        query = self.ncbi_query_edit.toPlainText()
        retmax = int(self.ncbi_retmax_edit.text())
        excel_file = f"{project_prefix}_geo_wide_metadata.xlsx"
        final_file = f"{project_prefix}_metadata_split_fixed_series_accession.xlsx"
        self.append_log("Starting NCBI data fetch...")
        self.thread = WorkerThread(
            run_ncbi_search, query, retmax, excel_file, final_file,
            progress_callback=self.update_progress
        )
        self.thread.finished_signal.connect(self.on_finished)
        self.thread.start()

    def handle_geo(self):
        project_prefix = self.project_name_edit.text()
        input_file = f"{project_prefix}_metadata_split_fixed_series_accession.xlsx"
        output_folder = f"{project_prefix}_outputs"
        self.append_log("Starting GEO data processing...")
        self.thread = WorkerThread(
            run_process_geo_data, input_file, output_folder,
            progress_callback=self.update_progress
        )
        self.thread.finished_signal.connect(self.on_finished)
        self.thread.start()

    def handle_batch(self):
        if not openai.api_key or not openai.api_key.strip():
            QMessageBox.warning(self, "Missing API Key", "Please set your OpenAI API Key in the settings.")
            self.open_settings_dialog()
            if not openai.api_key or not openai.api_key.strip():
                return
        global client
        client = openai
        project_prefix = self.project_name_edit.text()
        system_prompt = self.batch_system_prompt_edit.toPlainText()
        description = self.batch_user_prompt_edit.toPlainText().strip()
        questions = [pair[1].text().strip() for pair in self.batch_question_edits if pair[1].text().strip()]
        if questions:
            questions_text = "\n".join(f"{i+1}. {q}" for i, q in enumerate(questions))
            combined_user_prompt = f"{description}\n\nQuestions:\n{questions_text}"
        else:
            combined_user_prompt = description
        input_directory = f"{project_prefix}_outputs"
        batch_file_path = f"{project_prefix}_batch_requests.jsonl"
        output_prefix = f"{project_prefix}_batch_output"
        self.append_log("Starting batch job creation...")
        self.thread = WorkerThread(
            run_create_batch_job,
            system_prompt,
            combined_user_prompt,
            input_directory,
            batch_file_path,
            output_prefix,
            progress_callback=self.update_progress
        )
        self.thread.finished_signal.connect(self.on_finished)
        self.thread.start()

    def handle_parse(self):
        s = QSettings("MyCompany", "NCBIApp")
        last_batch_id = s.value("last_batch_id", "", type=str)
        output_prefix = s.value("last_output_prefix", "", type=str)
        if not last_batch_id:
            self.append_log("No previous batch job found.")
            return
        self.append_log(f"Retrieving status for batch job {last_batch_id}...")
        
        global client
        if client is None:
            api_key = s.value("openai_api_key", "", type=str)
            if not api_key.strip():
                QMessageBox.warning(self, "Missing API Key", "Please set your OpenAI API Key in the settings.")
                self.open_settings_dialog()
                api_key = s.value("openai_api_key", "", type=str)
            else: 
                openai.api_key = api_key
        client = openai
        response = client.batches.retrieve(last_batch_id)
        status = response.status
        self.append_log(f"Current batch job status: {status}")
        if status != "completed":
            self.append_log("Batch job is not yet completed. Please wait.")
            return
        if response.output_file_id:
            jsonl_file = f"{output_prefix}.jsonl"
            csv_output = f"{output_prefix}.csv"
            download_file(response.output_file_id, jsonl_file)
            self.append_log(f"Batch output downloaded to {jsonl_file}. Parsing output...")
            parse_batch_output(jsonl_file, csv_output)
            self.append_log(f"Batch output parsed and saved to {csv_output}.")
            self.reveal_file_in_explorer(csv_output)
        else:
            self.append_log("No output file found for the completed batch job.")

    def on_finished(self, message):
        self.append_log(message)
        QMessageBox.information(self, "Task Finished", message)
    
    def open_settings_dialog(self):
        dialog = SettingsDialog(self)
        if dialog.exec():
            self.append_log("Settings updated.")
    
    def add_question_field(self):
        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(20, 0, 0, 0)
        
        line_edit = QLineEdit()
        line_edit.setPlaceholderText("Enter question")
        remove_button = QPushButton("✖")
        remove_button.setFlat(True)
        remove_button.setStyleSheet("""
            QPushButton {
                background: transparent;
                border: none;
                font-size: 16px;
                color: gray;
                padding: 4px;
            }
            QPushButton:hover {
                color: red;
            }
        """)
        remove_button.setMaximumWidth(30)
        remove_button.clicked.connect(lambda: self.remove_question_field(container, line_edit))
        layout.addWidget(line_edit)
        layout.addWidget(remove_button)
        
        self.batch_questions_layout.addWidget(container)
        self.batch_question_edits.append((container, line_edit))
        
    def remove_question_field(self, container, line_edit):
        self.batch_questions_layout.removeWidget(container)
        container.deleteLater()
        self.batch_question_edits = [pair for pair in self.batch_question_edits if pair[1] != line_edit]
    
    def toggle_dark_mode(self):
        self.dark_mode = not self.dark_mode
        apply_styles(self.dark_mode)
        if self.dark_mode:
            self.dark_mode_button.setText("☀️")
        else:
            self.dark_mode_button.setText("🌙")

    def handle_project_name_changed(self, text):
        if " " in text:
            new_text = text.replace(" ", "-")
            self.project_name_edit.blockSignals(True)
            self.project_name_edit.setText(new_text)
            self.project_name_edit.blockSignals(False)
    
    def reveal_file_in_explorer(self, path):
        abs_path = os.path.abspath(path)
        if sys.platform == "darwin":
            subprocess.run(["open", "-R", abs_path])
        elif sys.platform == "win32":
            subprocess.run(["explorer", "/select,", abs_path])
        else:
            subprocess.run(["xdg-open", os.path.dirname(abs_path)])
    
    @Slot(str)
    def append_log(self, msg):
        self.log_text.append(msg)

# Appearance
def is_macos_dark_mode():
    if sys.platform != "darwin":
        return False
    try:
        result = subprocess.run(["defaults", "read", "-g", "AppleInterfaceStyle"],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout.strip() == "Dark"
    except Exception:
        return False

def apply_styles(dark_mode):
    if dark_mode:
        style = """
            * {
                font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
                font-size: 13px;
                color: #FFFFFF;
            }
            QMainWindow { background-color: #2B2B2B; }
            QPushButton {
                background: transparent;
                border: none;
                color: white;
                padding: 8px 16px;
                border-radius: 6px;
            }
            QPushButton:hover { color: #0A84FF; }
            QLineEdit, QTextEdit {
                background-color: #3C3C3C;
                border: 1px solid #555;
                padding: 6px;
                border-radius: 4px;
                color: #FFFFFF;
            }
            QLineEdit:focus, QTextEdit:focus { border: 1px solid #0A84FF; }
            QGroupBox {
                background-color: #2D2D2D;
                border: 1px solid #555;
                border-radius: 6px;
                margin-top: 20px;
                padding: 10px;
                font-weight: bold; 
            }
            QGroupBox:title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }
            QProgressBar {
                border: 1px solid #555;
                text-align: center;
                background-color: #3C3C3C;
                color: #FFFFFF;
                border-radius: 4px;
            }
            QProgressBar::chunk {
                background-color: #0A84FF;
                border-radius: 4px;
            }
            QMessageBox, QDialog {
                background-color: #2B2B2B;
                color: #FFFFFF;
            }
            /* NEW: Modern macOS-like scroll bar */
            QScrollBar:vertical {
                background: transparent;
                width: 8px;
                margin: 0px;
            }
            QScrollBar::handle:vertical {
                background: #555;
                border-radius: 4px;
                min-height: 20px;
            }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                height: 0px;
            }
            QScrollBar:horizontal {
                background: transparent;
                height: 8px;
                margin: 0px;
            }
            QScrollBar::handle:horizontal {
                background: #555;
                border-radius: 4px;
                min-width: 20px;
            }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                width: 0px;
            }
        """
    else:
        style = """
            * {
                font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
                font-size: 13px;
                color: #000000;
            }
            QMainWindow { background-color: #F2F2F2; }
            QPushButton {
                background: transparent;
                border: none;
                color: #007AFF;
                padding: 8px 16px;
                border-radius: 6px;
            }
            QPushButton:hover { color: #339CFF; }
            QLineEdit, QTextEdit {
                background-color: #FFFFFF;
                border: 1px solid #CCC;
                padding: 6px;
                border-radius: 4px;
                color: #000000;
            }
            QLineEdit:focus, QTextEdit:focus { border: 1px solid #007AFF; }
            QGroupBox {
                background-color: #EFEFEF;
                border: 1px solid #CCC;
                border-radius: 6px;
                margin-top: 20px;
                padding: 10px;
                font-weight: bold; 
            }
            QGroupBox:title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }
            QProgressBar {
                border: 1px solid #CCC;
                text-align: center;
                background-color: #FFFFFF;
                color: #000000;
                border-radius: 4px;
            }
            QProgressBar::chunk {
                background-color: #007AFF;
                border-radius: 4px;
            }
            QMessageBox, QDialog {
                background-color: #F2F2F2;
                color: #000000;
            }
            /* NEW: Modern macOS-like scroll bar */
            QScrollBar:vertical {
                background: transparent;
                width: 8px;
                margin: 0px;
            }
            QScrollBar::handle:vertical {
                background: #888;
                border-radius: 4px;
                min-height: 20px;
            }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                height: 0px;
            }
            QScrollBar:horizontal {
                background: transparent;
                height: 8px;
                margin: 0px;
            }
            QScrollBar::handle:horizontal {
                background: #888;
                border-radius: 4px;
                min-width: 20px;
            }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                width: 0px;
            }
        """
    QApplication.instance().setStyleSheet(style)

if __name__ == "__main__":
    app = QApplication([])
    app.setStyle("Fusion")
    dark_mode = is_macos_dark_mode() if sys.platform == "darwin" else False
    apply_styles(dark_mode)
    window = MainWindow()
    window.show()
    app.exec()