#!/usr/bin/env python3
import re
import requests
from bs4 import BeautifulSoup
import spacy

# Load SpaCy model
nlp = spacy.load("en_core_web_sm")

# TCGA Cancer Type Abbreviations with common synonyms
cancer_types = {
    "LAML": ["acute myeloid leukemia"],
    "ACC": ["adrenocortical carcinoma"],
    "BLCA": ["bladder urothelial carcinoma"],
    "BCC": ["basal cell carcinoma","squamous cell carcinoma"], # Skin
    "LGG": ["brain lower grade glioma"],
    "BRCA": ["breast invasive carcinoma"],
    "CESC": ["cervical squamous cell carcinoma", "endocervical adenocarcinoma"],
    "CHOL": ["cholangiocarcinoma"],
    "LCML": ["chronic myelogenous leukemia"],
    "COAD": ["colon adenocarcinoma", "colorectal cancer", "colorectal carcinoma"],
    "CNTL": ["controls"],
    "ESCA": ["esophageal carcinoma"],
    "FPPP": ["ffpe pilot phase ii"],
    "GBM": ["glioblastoma multiforme"],
    "HNSC": ["head and neck squamous cell carcinoma", "head and neck","nasopharyngeal"],
    "KICH": ["kidney chromophobe"],
    "KIRC": ["kidney renal clear cell carcinoma"],
    "KIRP": ["kidney renal papillary cell carcinoma"],
    "LIHC": ["liver hepatocellular carcinoma", "hepatocellular carcinoma"],
    "LUAD": ["lung adenocarcinoma"],
    "LUSC": ["lung squamous cell carcinoma"],
    "DLBC": ["lymphoid neoplasm diffuse large b-cell lymphoma"],
    "MESO": ["mesothelioma"],
    "MISC": ["miscellaneous"],
    "OV": ["ovarian serous cystadenocarcinoma"],
    "PAAD": ["pancreatic adenocarcinoma"],
    "PCPG": ["pheochromocytoma", "paraganglioma"],
    "PRAD": ["prostate adenocarcinoma"],
    "READ": ["rectum adenocarcinoma"],
    "SARC": ["sarcoma"],
    "SKCM": ["skin cutaneous melanoma"],
    "STAD": ["stomach adenocarcinoma", "gastric adenocarcinoma", "gastric tumor", "gastric"],
    "TGCT": ["testicular germ cell tumors"],
    "THYM": ["thymoma"],
    "THCA": ["thyroid carcinoma"],
    "UCS": ["uterine carcinosarcoma"],
    "UCEC": ["uterine corpus endometrial carcinoma"],
    "UVM": ["uveal melanoma"]
}

def get_cancer_type(text):
    text = text.lower()
    for acronym, names in cancer_types.items():
        for name in names:
            if name in text:
                # print(f"Matched cancer type '{name}' as '{acronym}' in text.")  # Debugging output
                return acronym
    #print(f"No match found for cancer type in text: {text}")  # Debugging output
    return "UNKNOWN"

def classify_unknown(text):
    text = text.lower()
    if "pancancer" in text:
        return "PANCANCER"
    if any(term in text for term in ["tumor", "tumour", "carcinoma", "cancer"]):
        return "CANCER"
    if any(term in text for term in ["inflammatory bowel disease", "rheumatoid arthritis"]):
        return "INFLAM"
    return "UNKNOWN"

def parse_pubmed(pubmed_id):
    # Prepend PubMed URL
    url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
    
    # Fetch the page content
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Extract title and year
    title = soup.find('h1', class_='heading-title').get_text(strip=True)
    citation_text = soup.find('span', class_='cit').get_text()
    year = re.search(r'\b\d{4}\b', citation_text).group(0)
    
    # Extract abstract
    abstract = soup.find('div', class_='abstract').get_text(strip=True)
    
    # Extract first author's last name
    first_author_name = soup.find('a', class_='full-name').get_text(strip=True)
    first_author_lastname = first_author_name.split(' ')[-1].upper()
    
    # Get the cancer type from title first, then abstract as fallback
    cancer_type = get_cancer_type(title)
    if cancer_type == "UNKNOWN":
        cancer_type = get_cancer_type(abstract)
    
    # Additional classification for UNKNOWN
    if cancer_type == "UNKNOWN":
        cancer_type = classify_unknown(abstract)
    
    # Combine into the desired format
    result = f"{cancer_type}_{year}_{pubmed_id}_{first_author_lastname}"
    return result

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <pubmed_id>")
        sys.exit(1)

    pubmed_id = sys.argv[1]
    result = parse_pubmed(pubmed_id)
    print(result)

