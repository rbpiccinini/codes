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

def getdata(campo,d0,d1,idx):
    # Connect to duckduckgo
    browser = mechanicalsoup.StatefulBrowser()
    browser.open("https://www.anp.gov.br/SITE/extras/consulta_petroleo_derivados/producao/consultaProdMensalHidrocarbonetos/default.asp", verify=False)
    
    # Fill-in the search form
    browser.select_form('form[action="default.asp"]')
    browser["txtDe"] = d0
    browser["txtAte"] = d1
    #browser['txtCampo'] = campo.encode("utf-8")
    response1 = browser.submit_selected()
    
    # Display the results
    page = browser.get_current_page()
    browser.select_form('form[action="exportar.asp"]')
    try:
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
    except:
        raise ValueError('Sem dados.')

#campos = pd.read_excel('producao.xlsx')['Campo'].tolist()
dfs=[]


#years = range(1950,2020,10)
#for campo in campos:
#    print(campo)
#    for year in years:
#        d0 ='01/'+str(year)
#        d1 ='12/'+str(year+9)
#        try:
#            df = getdata(campo, d0,d1, './Campos_Produção/'+campo+str(year))
#            print(campo + ' - '+ d0+' a '+d1+' - Linhas = '+ str(len(df)))
#            dfs.append(df)
#            time.sleep(0.5)
#        except:
#            print(campo + ' - '+ d0+' a '+d1+' - Sem dados.')
#            continue
        
years = range(1950,2020,1)
for year in years:
    d0 ='01/'+str(year)
    d1 ='12/'+str(year)
    try:
        df = getdata('', d0,d1, './Campos_Produção/'+str(year))
        print(' - '+ d0+' a '+d1+' - Linhas = '+ str(len(df)))
        dfs.append(df)
        time.sleep(0.5)
    except:
        print(' - '+ d0+' a '+d1+' - Sem dados.')
        continue
    
pocos=pd.concat(dfs)
pocos.to_csv('campos_prod.csv', sep=';')
