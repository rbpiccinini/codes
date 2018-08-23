import mechanicalsoup
import requests
import pandas as pd
import numpy as np
import time
import urllib3
import datetime

# disable ssl warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
# disable ssl verification
session = requests.Session()
session.verify = False

def getdata(d0,d1,idx):
    # Connect to duckduckgo
    browser = mechanicalsoup.StatefulBrowser()
    browser.open("https://www.anp.gov.br/SITE/extras/consulta_petroleo_derivados/producao/consultaProdCamposEmProducao/default.asp", verify=False)
    
    # Fill-in the search form
    browser.select_form('form[action="default.asp"]')
    browser["txtDe"] = d0
    browser["txtAte"] = d1
    response1 = browser.submit_selected()
    
    # Display the results
    page = browser.get_current_page()
    browser.select_form('form[action="exportar.asp"]')
    browser["Sim"] = 'Sim'
    response2 = browser.submit_selected()
    
    #response = browser.open('https://www.anp.gov.br/SITE/extras/consulta_petroleo_derivados/exploracao/consultaExploPocosPerfurados/exportar.asp')
    browser.select_form()
    browser["Sim"] = 'Sim'
    response3 = browser.submit_selected()
    f = open(idx+'.html','w')
    f.write(response3.text)
    f.close()
    
    df = pd.read_html(idx+'.html', match='.+', flavor=None, header=0, index_col=None, skiprows=None, attrs=None, parse_dates=True, tupleize_cols=None, thousands='.', encoding='utf8', decimal=',', converters=None, na_values=None, keep_default_na=True, displayed_only=True)
    return df[-1]


dfs=[]
years = range(2010,2019,1)

for year in years:
    d0 ='01/'+str(year)
    d1 ='12/'+str(year)
    df = getdata(d0,d1, './Campos_Produção/'+str(year))
    print(str(year) + ' - Linhas = '+ str(df.shape[0]))
    dfs.append(df)
    time.sleep(0.5)
    
pocos=pd.concat(dfs)
pocos.to_csv('campos_prod.csv', sep=';')
