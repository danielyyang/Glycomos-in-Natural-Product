from Bio import Entrez
import logging
import time
from typing import List

# 配置日志 (Configure logging)
logger = logging.getLogger(__name__)

class LiteratureResearcher:
    """
    负责文献挖掘与生物活性提取 (Responsible for literature mining and bioactivity extraction)
    """
    
    def __init__(self, email: str):
        # Entrez requires an email
        Entrez.email = email

    def search_pubmed(self, query: str, max_results: int = 10) -> List[str]:
        """
        在 PubMed 中搜索相关文献摘要 (Search PubMed for relevant abstracts)
        """
        abstracts = []
        try:
            logger.info(f"Searching PubMed for: {query}")
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if id_list:
                handle = Entrez.efetch(db="pubmed", id=id_list, rettype="abstract", retmode="text")
                text_data = handle.read()
                handle.close()
                abstracts.append(text_data)
                
            time.sleep(1) # API limit
            
        except Exception as e:
            logger.error(f"Error searching PubMed: {e}")
            
        return abstracts

    def extract_keywords(self, text: str) -> List[str]:
        """
        简单的关键词提取 (Simple keyword extraction)
        """
        keywords = ["anti-inflammatory", "antitumor", "surfactant", "cytotoxic"]
        found = [k for k in keywords if k in text.lower()]
        return found

if __name__ == "__main__":
    # Example email, user should replace
    researcher = LiteratureResearcher("example@email.com")
    results = researcher.search_pubmed("Solamargine bioactivity")
    for res in results:
        print(researcher.extract_keywords(res))
