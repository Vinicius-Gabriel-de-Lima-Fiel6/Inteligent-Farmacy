import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import requests
import difflib
from chempy import balance_stoichiometry
from scipy.interpolate import make_interp_spline
from supabase import create_client, Client

# --- 1. CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(page_title="BioPharm Ultra 2026", layout="wide", page_icon="🧪")

# --- 2. CONEXÃO COM BANCO DE DADOS (SUPABASE) ---
try:
    url: str = st.secrets["SUPABASE_URL"]
    key: str = st.secrets["SUPABASE_KEY"]
    supabase: Client = create_client(url, key)
except Exception:
    st.error("⚠️ Erro: Configure as chaves do Supabase nas Secrets para ativar o banco em nuvem.")

# --- 3. FUNÇÕES DE API EXTERNA (PUBCHEM) ---
def busca_api_pubchem(termo):
    """Busca dados técnicos em tempo real via API do PubChem"""
    try:
        url_api = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{termo}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        res = requests.get(url_api, timeout=5)
        if res.status_code == 200:
            d = res.json()['PropertyTable']['Properties'][0]
            return {
                "formula": d.get('MolecularFormula'),
                "massa": d.get('MolecularWeight'),
                "nome": d.get('IUPACName')
            }
    except:
        return None

# --- 4. ESTILIZAÇÃO CSS CUSTOMIZADA ---
st.markdown("""
    <style>
    .main-header {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        padding: 2rem; border-radius: 15px; color: white; text-align: center; margin-bottom: 2rem;
    }
    .element-card {
        border-radius: 10px; padding: 15px; text-align: center; color: white;
        font-weight: bold; border: 1px solid rgba(255,255,255,0.1); height: 110px;
        box-shadow: 2px 2px 10px rgba(0,0,0,0.1);
    }
    </style>
    """, unsafe_allow_html=True)

# --- 5. DADOS ESTÁTICOS (TABELA PERIÓDICA) ---
ELEMENTOS = {
    "H": {"Nome": "Hidrogênio", "Número": 1, "Massa": 1.008, "Grupo": "Não metais", "Período": "1", "Descrição": "Elemento mais leve e abundante."},
"He": {"Nome": "Hélio", "Número": 2, "Massa": 4.0026, "Grupo": "Gases nobres", "Período": "1", "Descrição": "Gás nobre, usado em balões."},
"Li": {"Nome": "Lítio", "Número": 3, "Massa": 6.94, "Grupo": "Metais alcalinos", "Período": "2", "Descrição": "Metal pouco reativo."},
"Be": {"Nome": "Berílio", "Número": 4, "Massa": 9.0122, "Grupo": "Metais alcalino-terrosos", "Período": "2", "Descrição": "Carga para estabilidade é +2."},
"B": {"Nome": "Boro", "Número": 5, "Massa": 10.81, "Grupo": "Semimetais", "Período": "2", "Descrição": "Precisa perder 3 elétrons para ficar estável."},
"C": {"Nome": "Carbono", "Número": 6, "Massa": 12.011, "Grupo": "Não metais", "Período": "2", "Descrição": "Tetravalente e forma cadeias carbônicas."},
"N": {"Nome": "Nitrogênio", "Número": 7, "Massa": 14.007, "Grupo": "Não metais", "Período": "2", "Descrição": "Elemento fundamental para a vida."},
"O": {"Nome": "Oxigênio", "Número": 8, "Massa": 15.999, "Grupo": "Não metais", "Período": "2", "Descrição": "Essencial para a respiração."},
"F": {"Nome": "Flúor", "Número": 9, "Massa": 18.998, "Grupo": "Halogênios", "Período": "2", "Descrição": "Mais eletronegativo da tabela."},
"Ne": {"Nome": "Neônio", "Número": 10, "Massa": 20.180, "Grupo": "Gases nobres", "Período": "2", "Descrição": "Usado em neons e luzes."},
"Na": {"Nome": "Sódio", "Número": 11, "Massa": 22.990, "Grupo": "Metais alcalinos", "Período": "3", "Descrição": "Muito eletropositivo e seus eletrólitos conduzem eletricidade."},
"Mg": {"Nome": "Magnésio", "Número": 12, "Massa": 24.305, "Grupo": "Metais alcalino-terrosos", "Período": "3", "Descrição": "Pode formar o leite de magnésia."},
"Al": {"Nome": "Alumínio", "Número": 13, "Massa": 26.982, "Grupo": "Outros metais", "Período": "3", "Descrição": "Metal forte e muito maleável."},
"Si": {"Nome": "Silício", "Número": 14, "Massa": 28.085, "Grupo": "Semimetais", "Período": "3", "Descrição": "Usado em semicondutores."},
"P": {"Nome": "Fósforo", "Número": 15, "Massa": 30.974, "Grupo": "Não metais", "Período": "3", "Descrição": "Fundamental à vida e precisa de 3 elétrons para a estabilidade."},
"S": {"Nome": "Enxofre", "Número": 16, "Massa": 32.06, "Grupo": "Não metais", "Período": "3", "Descrição": "Fundamental à vida e precisa de 2 elétrons."},
"Cl": {"Nome": "Cloro", "Número": 17, "Massa": 35.45, "Grupo": "Halogênios", "Período": "3", "Descrição": "Altamente eletronegativo e precisa de 1 elétron."},
"Ar": {"Nome": "Argônio", "Número": 18, "Massa": 39.948, "Grupo": "Gases nobres", "Período": "3", "Descrição": "Gás nobre, usado em lâmpadas e fotografia."},
"K": {"Nome": "Potássio", "Número": 19, "Massa": 39.098, "Grupo": "Metais alcalinos", "Período": "4", "Descrição": "Metal alcalino altamente reativo, essencial para funções celulares."},
"Ca": {"Nome": "Cálcio", "Número": 20, "Massa": 40.078, "Grupo": "Metais alcalino-terrosos", "Período": "4", "Descrição": "Importante para ossos e dentes."},
"Sc": {"Nome": "Escândio", "Número": 21, "Massa": 44.956, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Usado em ligas leves."},
"Ti": {"Nome": "Titânio", "Número": 22, "Massa": 47.867, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Forte e resistente à corrosão."},
"V": {"Nome": "Vanádio", "Número": 23, "Massa": 50.942, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Fortalece aço."},
"Cr": {"Nome": "Cromo", "Número": 24, "Massa": 51.996, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Usado em cromagem e ligas."},
"Mn": {"Nome": "Manganês", "Número": 25, "Massa": 54.938, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Importante para ligas de aço."},
"Fe": {"Nome": "Ferro", "Número": 26, "Massa": 55.845, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Essencial na hemoglobina."},
"Co": {"Nome": "Cobalto", "Número": 27, "Massa": 58.933, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Usado em ímãs e baterias."},
"Ni": {"Nome": "Níquel", "Número": 28, "Massa": 58.693, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Usado em ligas e moedas."},
"Cu": {"Nome": "Cobre", "Número": 29, "Massa": 63.546, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Excelente condutor elétrico."},
"Zn": {"Nome": "Zinco", "Número": 30, "Massa": 65.38, "Grupo": "Metais de transição", "Período": "4", "Descrição": "Galvanização e essencial ao organismo."},
"Ga": {"Nome": "Gálio", "Número": 31, "Massa": 69.723, "Grupo": "Outros metais", "Período": "4", "Descrição": "Derrete na mão, usado em eletrônicos."},
"Ge": {"Nome": "Germânio", "Número": 32, "Massa": 72.63, "Grupo": "Semimetais", "Período": "4", "Descrição": "Usado em semicondutores."},
"As": {"Nome": "Arsênio", "Número": 33, "Massa": 74.922, "Grupo": "Semimetais", "Período": "4", "Descrição": "Tóxico, usado em pesticidas."},
"Se": {"Nome": "Selênio", "Número": 34, "Massa": 78.971, "Grupo": "Não metais", "Período": "4", "Descrição": "Essencial em pequenas quantidades."},
"Br": {"Nome": "Bromo", "Número": 35, "Massa": 79.904, "Grupo": "Halogênios", "Período": "4", "Descrição": "Líquido, usado em retardadores de chama."},
"Kr": {"Nome": "Criptônio", "Número": 36, "Massa": 83.798, "Grupo": "Gases nobres", "Período": "4", "Descrição": "Usado em lâmpadas e fotografia."},
"Rb": {"Nome": "Rubídio", "Número": 37, "Massa": 85.468, "Grupo": "Metais alcalinos", "Período": "5", "Descrição": "Altamente reativo, usado em pesquisas."},
"Sr": {"Nome": "Estrôncio", "Número": 38, "Massa": 87.62, "Grupo": "Metais alcalino-terrosos", "Período": "5", "Descrição": "Fogos de artifício e ligas metálicas."},
"Y": {"Nome": "Ítrio", "Número": 39, "Massa": 88.906, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Usado em LEDs e supercondutores."},
"Zr": {"Nome": "Zircônio", "Número": 40, "Massa": 91.224, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Resistente à corrosão, usado em reatores."},
"Nb": {"Nome": "Nióbio", "Número": 41, "Massa": 92.906, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Usado para fortalecer aço e em supercondutores."},
"Mo": {"Nome": "Molibdênio", "Número": 42, "Massa": 95.95, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Essencial em ligas e enzimas."},
"Tc": {"Nome": "Tecnécio", "Número": 43, "Massa": 98, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Radioativo, usado em medicina nuclear."},
"Ru": {"Nome": "Rutênio", "Número": 44, "Massa": 101.07, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Catalisador e ligas elétricas."},
"Rh": {"Nome": "Ródio", "Número": 45, "Massa": 102.91, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Catalisadores automotivos."},
"Pd": {"Nome": "Paládio", "Número": 46, "Massa": 106.42, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Joalheria e catalisadores."},
"Ag": {"Nome": "Prata", "Número": 47, "Massa": 107.87, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Melhor condutor elétrico."},
"Cd": {"Nome": "Cádmio", "Número": 48, "Massa": 112.41, "Grupo": "Metais de transição", "Período": "5", "Descrição": "Baterias e revestimentos."},
"In": {"Nome": "Índio", "Número": 49, "Massa": 114.82, "Grupo": "Outros metais", "Período": "5", "Descrição": "Telas sensíveis ao toque."},
"Sn": {"Nome": "Estanho", "Número": 50, "Massa": 118.71, "Grupo": "Outros metais", "Período": "5", "Descrição": "Bronze, soldas."},
"Sb": {"Nome": "Antimônio", "Número": 51, "Massa": 121.76, "Grupo": "Semimetais", "Período": "5", "Descrição": "Retardadores de chama e ligas."},
"Te": {"Nome": "Telúrio", "Número": 52, "Massa": 127.60, "Grupo": "Semimetais", "Período": "5", "Descrição": "Ligas metálicas e semicondutores."},
"I": {"Nome": "Iodo", "Número": 53, "Massa": 126.90, "Grupo": "Halogênios", "Período": "5", "Descrição": "Função da tireoide, antissépticos."},
"Xe": {"Nome": "Xenônio", "Número": 54, "Massa": 131.29, "Grupo": "Gases nobres", "Período": "5", "Descrição": "Lâmpadas e anestesia."},
"Cs": {"Nome": "Césio", "Número": 55, "Massa": 132.91, "Grupo": "Metais alcalinos", "Período": "6", "Descrição": "Relógios atômicos."},
"Ba": {"Nome": "Bário", "Número": 56, "Massa": 137.33, "Grupo": "Metais alcalino-terrosos", "Período": "6", "Descrição": "Radiologia, fogos de artifício."},
"La": {"Nome": "Lantânio", "Número": 57, "Massa": 138.91, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Inicia os lantanídeos, lentes ópticas."},
"Ce": {"Nome": "Cério", "Número": 58, "Massa": 140.12, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Catalisadores, polidores."},
"Pr": {"Nome": "Praseodímio", "Número": 59, "Massa": 140.91, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Ímãs e ligas aeronáuticas."},
"Nd": {"Nome": "Neodímio", "Número": 60, "Massa": 144.24, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Ímãs potentes."},
"Pm": {"Nome": "Promécio", "Número": 61, "Massa": 145, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Radioativo, baterias nucleares."},
"Sm": {"Nome": "Samário", "Número": 62, "Massa": 150.36, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Ímãs e lasers."},
"Eu": {"Nome": "Európio", "Número": 63, "Massa": 151.96, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Fósforos de telas e lâmpadas."},
"Gd": {"Nome": "Gadolínio", "Número": 64, "Massa": 157.25, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Ressonância magnética."},
"Tb": {"Nome": "Térbio", "Número": 65, "Massa": 158.93, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Fósforos verdes e eletrônicos."},
"Dy": {"Nome": "Disprósio", "Número": 66, "Massa": 162.50, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Ímãs e lasers."},
"Ho": {"Nome": "Hólmio", "Número": 67, "Massa": 164.93, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Ímãs e aplicações nucleares."},
"Er": {"Nome": "Érbio", "Número": 68, "Massa": 167.26, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Fibras ópticas e lasers médicos."},
"Tm": {"Nome": "Túlio", "Número": 69, "Massa": 168.93, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Lasers portáteis."},
"Yb": {"Nome": "Itérbio", "Número": 70, "Massa": 173.05, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Relógios atômicos e materiais especiais."},
"Lu": {"Nome": "Lutécio", "Número": 71, "Massa": 174.97, "Grupo": "Lantanídeos", "Período": "6", "Descrição": "Tomografia e catálise."},
"Hf": {"Nome": "Háfnio", "Número": 72, "Massa": 178.49, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Ligas de alta temperatura."},
"Ta": {"Nome": "Tântalo", "Número": 73, "Massa": 180.95, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Eletrônicos e implantes."},
"W": {"Nome": "Tungstênio", "Número": 74, "Massa": 183.84, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Mais alto ponto de fusão."},
"Re": {"Nome": "Rênio", "Número": 75, "Massa": 186.21, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Ligas e catalisadores."},
"Os": {"Nome": "Ósmio", "Número": 76, "Massa": 190.23, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Metal mais denso."},
"Ir": {"Nome": "Irídio", "Número": 77, "Massa": 192.22, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Equipamentos médicos, resistente à corrosão."},
"Pt": {"Nome": "Platina", "Número": 78, "Massa": 195.08, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Joias e catalisadores."},
"Au": {"Nome": "Ouro", "Número": 79, "Massa": 196.97, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Metal precioso e maleável."},
"Hg": {"Nome": "Mercúrio", "Número": 80, "Massa": 200.59, "Grupo": "Metais de transição", "Período": "6", "Descrição": "Único metal líquido à temperatura ambiente."},
"Tl": {"Nome": "Tálio", "Número": 81, "Massa": 204.38, "Grupo": "Outros metais", "Período": "6", "Descrição": "Tóxico, usado em eletrônicos."},
"Pb": {"Nome": "Chumbo", "Número": 82, "Massa": 207.2, "Grupo": "Outros metais", "Período": "6", "Descrição": "Denso, usado em proteção contra radiação."},
"Bi": {"Nome": "Bismuto", "Número": 83, "Massa": 208.98, "Grupo": "Outros metais", "Período": "6", "Descrição": "Menos tóxico que o chumbo."},
"Po": {"Nome": "Polônio", "Número": 84, "Massa": 209, "Grupo": "Semimetais", "Período": "6", "Descrição": "Radioativo, fontes de calor."},
"At": {"Nome": "Astato", "Número": 85, "Massa": 210, "Grupo": "Halogênios", "Período": "6", "Descrição": "Raro e radioativo."},
"Rn": {"Nome": "Radônio", "Número": 86, "Massa": 222, "Grupo": "Gases nobres", "Período": "6", "Descrição": "Radioativo, perigoso em ambientes fechados."},
"Fr": {"Nome": "Frâncio", "Número": 87, "Massa": 223, "Grupo": "Metais alcalinos", "Período": "7", "Descrição": "Extremamente raro e radioativo."},
"Ra": {"Nome": "Rádio", "Número": 88, "Massa": 226, "Grupo": "Metais alcalino-terrosos", "Período": "7", "Descrição": "Radioativo, usado em luminância antiga."},
"Ac": {"Nome": "Actínio", "Número": 89, "Massa": 227, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Inicia os actinídeos, altamente radioativo."},
"Th": {"Nome": "Tório", "Número": 90, "Massa": 232.04, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Potencial combustível nuclear."},
"Pa": {"Nome": "Protactínio", "Número": 91, "Massa": 231.04, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Radioativo, usado em pesquisas nucleares."},
"Np": {"Nome": "Neptúnio", "Número": 93, "Massa": 237, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Produzido em reatores nucleares, radioativo."},
"Pu": {"Nome": "Plutônio", "Número": 94, "Massa": 244, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Usado em armas nucleares e reatores."},
"Am": {"Nome": "Amerício", "Número": 95, "Massa": 243, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Detectores de fumaça."},
"Cm": {"Nome": "Cúrio", "Número": 96, "Massa": 247, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Fonte de energia espacial."},
"Bk": {"Nome": "Berquélio", "Número": 97, "Massa": 247, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Usado em pesquisa nuclear."},
"Cf": {"Nome": "Califórnio", "Número": 98, "Massa": 251, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Fonte de nêutrons."},
"Es": {"Nome": "Einstênio", "Número": 99, "Massa": 252, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Produzido em explosões nucleares."},
"Fm": {"Nome": "Férmio", "Número": 100, "Massa": 257, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Usado em estudos científicos."},
"Md": {"Nome": "Mendelévio", "Número": 101, "Massa": 258, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Pesquisa química nuclear."},
"No": {"Nome": "Nobélio", "Número": 102, "Massa": 259, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Experimentos laboratoriais."},
"Lr": {"Nome": "Laurêncio", "Número": 103, "Massa": 266, "Grupo": "Actinídeos", "Período": "7", "Descrição": "Elemento sintético radioativo."},
"Rf": {"Nome": "Rutherfórdio", "Número": 104, "Massa": 267, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Elemento sintético para pesquisa."},
"Db": {"Nome": "Dúbnio", "Número": 105, "Massa": 268, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Radioativo, instável."},
"Sg": {"Nome": "Seabórgio", "Número": 106, "Massa": 271, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Usado apenas em pesquisa."},
"Bh": {"Nome": "Bóhrio", "Número": 107, "Massa": 270, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Instável e sintético."},
"Hs": {"Nome": "Hássio", "Número": 108, "Massa": 277, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Pesquisado em laboratórios nucleares."},
"Mt": {"Nome": "Meitnério", "Número": 109, "Massa": 278, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Superpesado e sintético."},
"Ds": {"Nome": "Darmstádio", "Número": 110, "Massa": 281, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Meia-vida muito curta."},
"Rg": {"Nome": "Roentgênio", "Número": 111, "Massa": 282, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Elemento radioativo sintético."},
"Cn": {"Nome": "Copernício", "Número": 112, "Massa": 285, "Grupo": "Metais de transição", "Período": "7", "Descrição": "Altamente instável."},
"Nh": {"Nome": "Nihônio", "Número": 113, "Massa": 286, "Grupo": "Outros metais", "Período": "7", "Descrição": "Elemento sintético."},
"Fl": {"Nome": "Fleróvio", "Número": 114, "Massa": 289, "Grupo": "Outros metais", "Período": "7", "Descrição": "Superpesado, sintético."},
"Mc": {"Nome": "Moscóvio", "Número": 115, "Massa": 290, "Grupo": "Outros metais", "Período": "7", "Descrição": "Meia-vida curta."},
"Lv": {"Nome": "Livermório", "Número": 116, "Massa": 293, "Grupo": "Outros metais", "Período": "7", "Descrição": "Elemento instável."},
"Ts": {"Nome": "Tenessino", "Número": 117, "Massa": 294, "Grupo": "Halogênios", "Período": "7", "Descrição": "Superpesado e sintético."},
"Og": {"Nome": "Oganessônio", "Número": 118, "Massa": 294, "Grupo": "Gases nobres", "Período": "7", "Descrição": "Altamente radioativo e sintético."},
"Uue": {"Nome": "Ununennium", "Número": 119, "Massa": 315, "Grupo": "Metais alcalinos", "Período": "8", "Descrição": "Previsto como metal alcalino superpesado."},
"Ubn": {"Nome": "Unbinilium", "Número": 120, "Massa": 320, "Grupo": "Metais alcalino-terrosos", "Período": "8", "Descrição": "Previsto como metal alcalino-terroso."},
"Ubu": {"Nome": "Unbiunium", "Número": 121, "Massa": 326, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Primeiro dos superactinídeos."},
"Ubb": {"Nome": "Unbibium", "Número": 122, "Massa": 328, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Hipotético do grupo 4."},
"Ubt": {"Nome": "Unbitrium", "Número": 123, "Massa": 330, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Propriedades desconhecidas."},
"Ubq": {"Nome": "Unbiquadium", "Número": 124, "Massa": 332, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Ainda não sintetizado."},
"Ubp": {"Nome": "Unbipentium", "Número": 125, "Massa": 334, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Potencial de propriedades únicas."},
"Ubh": {"Nome": "Unbihexium", "Número": 126, "Massa": 336, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Previsto como altamente estável."},
"Ubs": {"Nome": "Unbiseptium", "Número": 127, "Massa": 338, "Grupo": "Elementos superpesados", "Período": "8", "Descrição": "Totalmente teórico, sem dados experimentais."},
 "U":{"Nome":"Urânio","Número":92,"Massa":238.0289,"Grupo":"Actinídeos","Período":"7","Descrição":"Combustível nuclear, ogivas,blindagem,corantes de virdro e cerâmica."}

    
}

# --- 6. HEADER PRINCIPAL ---
st.markdown('<div class="main-header"><h1>BioPharm Ultra 2026</h1><p>Sistema Unificado: Banco de Dados, APIs e Cálculos Avançados</p></div>', unsafe_allow_html=True)

# --- 7. NAVEGAÇÃO POR ABAS ---
tabs = st.tabs(["💬 Chatbot Híbrido", "💎 Tabelas Químicas", "⚖️ Estequiometria & 3D", "📈 Gráficos de Solubilidade", "⚙️ Admin (Upload)"])

# --- ABA 1: CHATBOT (SUPABASE + PUBCHEM) ---
with tabs[0]:
    st.subheader("Assistente Virtual Inteligente")
    if "messages" not in st.session_state: st.session_state.messages = []
    for m in st.session_state.messages:
        with st.chat_message(m["role"]): st.write(m["content"])

    if prompt := st.chat_input("Perqunte sobre um composto ou conceito..."):
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"): st.write(prompt)
        
        resposta = "Não encontrei informações sobre isso."
        
        # 1ª Tentativa: Banco de Dados Próprio (Supabase)
        try:
            res_db = supabase.table("conhecimento").select("pergunta, resposta").execute()
            dados_locais = {item['pergunta']: item['resposta'] for item in res_db.data}
            matches = difflib.get_close_matches(prompt.lower(), dados_locais.keys(), n=1, cutoff=0.6)
            if matches:
                resposta = dados_locais[matches[0]]
            else:
                # 2ª Tentativa: API Externa (PubChem)
                dados_api = busca_api_pubchem(prompt)
                if dados_api:
                    resposta = f"🔍 **Resultado via API PubChem:**\n\n**Nome:** {dados_api['nome']}\n\n**Fórmula:** {dados_api['formula']}\n\n**Massa Molar:** {dados_api['massa']} g/mol"
        except:
            resposta = "Erro ao conectar aos serviços de dados."

        st.session_state.messages.append({"role": "assistant", "content": resposta})
        with st.chat_message("assistant"): st.write(resposta)

# --- ABA 2: TABELAS ---
with tabs[1]:
    modo_tab = st.radio("Selecione:", ["Tabela Periódica", "Kps (Solubilidade)"], horizontal=True)
    if modo_tab == "Tabela Periódica":
        cols = st.columns(4)
        for i, (simb, info) in enumerate(ELEMENTOS.items()):
            with cols[i % 4]:
                st.markdown(f'<div class="element-card" style="background:{info["cor"]}">{info["n"]}<br><span style="font-size:24px">{simb}</span><br><small>{info["m"]}</small></div>', unsafe_allow_html=True)
                st.button(f"Detalhes {simb}", key=f"btn_{simb}", on_click=lambda s=simb: st.toast(f"Categoria: {ELEMENTOS[s]['cat']}"))
    else:
        df_kps = pd.DataFrame([
    ["Brometo de chumbo(II)", "PbBr₂", "7,9 x 10⁻⁵"],
    ["Brometo de cobre(I)", "CuBr", "4,2 x 10⁻⁸"],
    ["Brometo de prata", "AgBr", "7,7 x 10⁻¹³"],
    ["Cloreto de prata", "AgCl", "1,6 x 10⁻¹⁰"],
    ["Hidróxido de alumínio", "Al(OH)₃", "1,1 x 10⁻³³"],
    ["Hidróxido de ferro(III)", "Fe(OH)₃", "4 x 10⁻³⁸"],
    ["Hidróxido de magnésio", "Mg(OH)₂", "1,8 x 10⁻¹¹"],
    ["Carbonato de cálcio", "CaCO₃", "3,36 x 10⁻⁹"],
    ["Sulfato de bário", "BaSO₄", "1,1 x 10⁻¹⁰"]
        ], columns=["Nome","Fórmula", "Kps"])
        st.table(df_kps)

# --- ABA 3: ESTEQUIOMETRIA & 3D ---
with tabs[2]:
    st.subheader("Balanceamento de Reações")
    reacao_input = st.text_input("Insira a reação (Ex: H2 + O2 -> H2O)")
    if st.button("Executar Balanço"):
        try:
            reag, prod = reacao_input.split("->")
            r_list = [x.strip() for x in reag.split("+")]
            p_list = [x.strip() for x in prod.split("+")]
            reac_bal, prod_bal = balance_stoichiometry(r_list, p_list)
            st.success(f"Equação Balanceada: {dict(reac_bal)} -> {dict(prod_bal)}")
        except: st.error("Erro na sintaxe da equação. Use 'A + B -> C'.")
    
    st.divider()
    st.subheader("Visualizador Molecular 3D (API)")
    comp_3d = st.text_input("Nome do composto para visualização:")
    if comp_3d:
        st.link_button(f"Abrir {comp_3d} no PubChem 3D", f"https://pubchem.ncbi.nlm.nih.gov/#query={comp_3d}")

# --- ABA 4: GRÁFICOS (PLOTLY) ---
with tabs[3]:
    st.subheader("Curvas de Solubilidade Interativas")
    c1, c2 = st.columns([1, 2])
    with c1:
        comp_nome = st.text_input("Nome do Sal", "KNO3")
        temp_vals = st.text_input("Temps (°C)", "0, 20, 40, 60, 80")
        sol_vals = st.text_input("Solubilidade (g/100g)", "13, 32, 64, 110, 169")
    
    try:
        x = np.array([float(i) for i in temp_vals.split(",")])
        y = np.array([float(i) for i in sol_vals.split(",")])
        x_smooth = np.linspace(x.min(), x.max(), 300)
        y_smooth = make_interp_spline(x, y, k=2)(x_smooth)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x_smooth, y=y_smooth, name=comp_nome, line=dict(color='#10b981', width=4)))
        fig.add_trace(go.Scatter(x=x, y=y, mode='markers', name="Pontos Reais", marker=dict(size=10, color='white', line=dict(width=2, color='#10b981'))))
        fig.update_layout(template="dark", xaxis_title="Temperatura (°C)", yaxis_title="g/100g H₂O")
        c2.plotly_chart(fig, use_container_width=True)
    except: st.info("Insira valores numéricos separados por vírgula para gerar o gráfico.")

# --- ABA 5: ADMIN (UPLOAD SUPABASE) ---
with tabs[4]:
    st.subheader("Alimentar Inteligência Coletiva")
    st.info("O que você cadastrar aqui será salvo na nuvem e o Chatbot aprenderá imediatamente.")
    with st.form("admin_form"):
        pergunta_n = st.text_input("Pergunta ou Conceito:")
        resposta_n = st.text_area("Resposta Detalhada:")
        if st.form_submit_button("Fazer Upload para Nuvem"):
            if pergunta_n and resposta_n:
                try:
                    supabase.table("conhecimento").insert({"pergunta": pergunta_n.lower(), "resposta": resposta_n}).execute()
                    st.success("✅ Conhecimento integrado com sucesso!")
                except Exception as e: st.error(f"Erro ao salvar: {e}")
            else: st.warning("Preencha todos os campos.")
