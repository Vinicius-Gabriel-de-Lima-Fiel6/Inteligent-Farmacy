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
    "H": {"n": 1, "m": 1.008, "cat": "Não-metal", "cor": "#3b82f6"},
    "He": {"n": 2, "m": 4.002, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Li": {"n": 3, "m": 6.94, "cat": "Alcalino", "cor": "#f59e0b"},
    "Be": {"n": 4, "m": 9.012, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "C": {"n": 6, "m": 12.01, "cat": "Não-metal", "cor": "#3b82f6"},
    "O": {"n": 8, "m": 15.99, "cat": "Não-metal", "cor": "#3b82f6"},
    "Na": {"n": 11, "m": 22.98, "cat": "Alcalino", "cor": "#f59e0b"},
    "Fe": {"n": 26, "m": 55.84, "cat": "Transição", "cor": "#ef4444"},
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
        df_kps = pd.DataFrame([["AgCl", "1,6 x 10⁻¹⁰"], ["BaSO₄", "1,1 x 10⁻¹⁰"], ["CaCO₃", "3,36 x 10⁻⁹"]], columns=["Fórmula", "Kps"])
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
