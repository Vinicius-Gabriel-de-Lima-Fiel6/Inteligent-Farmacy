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
ELEMENTOS ={
        "H": {"n": 1, "m": 1.008, "cat": "Não-metal", "cor": "#3b82f6"},
    "He": {"n": 2, "m": 4.002, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Li": {"n": 3, "m": 6.94, "cat": "Alcalino", "cor": "#f59e0b"},
    "Be": {"n": 4, "m": 9.012, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "B": {"n": 5, "m": 10.81, "cat": "Semimetal", "cor": "#06b6d4"},
    "C": {"n": 6, "m": 12.01, "cat": "Não-metal", "cor": "#3b82f6"},
    "N": {"n": 7, "m": 14.007, "cat": "Não-metal", "cor": "#3b82f6"},
    "O": {"n": 8, "m": 15.99, "cat": "Não-metal", "cor": "#3b82f6"},
    "F": {"n": 9, "m": 18.998, "cat": "Halogênio", "cor": "#f43f5e"},
    "Ne": {"n": 10, "m": 20.18, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Na": {"n": 11, "m": 22.99, "cat": "Alcalino", "cor": "#f59e0b"},
    "Mg": {"n": 12, "m": 24.305, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "Al": {"n": 13, "m": 26.982, "cat": "Outro metal", "cor": "#6b7280"},
    "Si": {"n": 14, "m": 28.085, "cat": "Semimetal", "cor": "#06b6d4"},
    "P": {"n": 15, "m": 30.974, "cat": "Não-metal", "cor": "#3b82f6"},
    "S": {"n": 16, "m": 32.06, "cat": "Não-metal", "cor": "#3b82f6"},
    "Cl": {"n": 17, "m": 35.45, "cat": "Halogênio", "cor": "#f43f5e"},
    "Ar": {"n": 18, "m": 39.948, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "K": {"n": 19, "m": 39.098, "cat": "Alcalino", "cor": "#f59e0b"},
    "Ca": {"n": 20, "m": 40.078, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "Sc": {"n": 21, "m": 44.956, "cat": "Transição", "cor": "#ef4444"},
    "Ti": {"n": 22, "m": 47.867, "cat": "Transição", "cor": "#ef4444"},
    "V": {"n": 23, "m": 50.942, "cat": "Transição", "cor": "#ef4444"},
    "Cr": {"n": 24, "m": 51.996, "cat": "Transição", "cor": "#ef4444"},
    "Mn": {"n": 25, "m": 54.938, "cat": "Transição", "cor": "#ef4444"},
    "Fe": {"n": 26, "m": 55.845, "cat": "Transição", "cor": "#ef4444"},
    "Co": {"n": 27, "m": 58.933, "cat": "Transição", "cor": "#ef4444"},
    "Ni": {"n": 28, "m": 58.693, "cat": "Transição", "cor": "#ef4444"},
    "Cu": {"n": 29, "m": 63.546, "cat": "Transição", "cor": "#ef4444"},
    "Zn": {"n": 30, "m": 65.38, "cat": "Transição", "cor": "#ef4444"},
    "Ga": {"n": 31, "m": 69.723, "cat": "Outro metal", "cor": "#6b7280"},
    "Ge": {"n": 32, "m": 72.63, "cat": "Semimetal", "cor": "#06b6d4"},
    "As": {"n": 33, "m": 74.922, "cat": "Semimetal", "cor": "#06b6d4"},
    "Se": {"n": 34, "m": 78.971, "cat": "Não-metal", "cor": "#3b82f6"},
    "Br": {"n": 35, "m": 79.904, "cat": "Halogênio", "cor": "#f43f5e"},
    "Kr": {"n": 36, "m": 83.798, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Rb": {"n": 37, "m": 85.468, "cat": "Alcalino", "cor": "#f59e0b"},
    "Sr": {"n": 38, "m": 87.62, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "Y": {"n": 39, "m": 88.906, "cat": "Transição", "cor": "#ef4444"},
    "Zr": {"n": 40, "m": 91.224, "cat": "Transição", "cor": "#ef4444"},
    "Nb": {"n": 41, "m": 92.906, "cat": "Transição", "cor": "#ef4444"},
    "Mo": {"n": 42, "m": 95.95, "cat": "Transição", "cor": "#ef4444"},
    "Tc": {"n": 43, "m": 98.0, "cat": "Transição", "cor": "#ef4444"},
    "Ru": {"n": 44, "m": 101.07, "cat": "Transição", "cor": "#ef4444"},
    "Rh": {"n": 45, "m": 102.91, "cat": "Transição", "cor": "#ef4444"},
    "Pd": {"n": 46, "m": 106.42, "cat": "Transição", "cor": "#ef4444"},
    "Ag": {"n": 47, "m": 107.87, "cat": "Transição", "cor": "#ef4444"},
    "Cd": {"n": 48, "m": 112.41, "cat": "Transição", "cor": "#ef4444"},
    "In": {"n": 49, "m": 114.82, "cat": "Outro metal", "cor": "#6b7280"},
    "Sn": {"n": 50, "m": 118.71, "cat": "Outro metal", "cor": "#6b7280"},
    "Sb": {"n": 51, "m": 121.76, "cat": "Semimetal", "cor": "#06b6d4"},
    "Te": {"n": 52, "m": 127.6, "cat": "Semimetal", "cor": "#06b6d4"},
    "I": {"n": 53, "m": 126.9, "cat": "Halogênio", "cor": "#f43f5e"},
    "Xe": {"n": 54, "m": 131.29, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Cs": {"n": 55, "m": 132.91, "cat": "Alcalino", "cor": "#f59e0b"},
    "Ba": {"n": 56, "m": 137.33, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "La": {"n": 57, "m": 138.91, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Ce": {"n": 58, "m": 140.12, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Pr": {"n": 59, "m": 140.91, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Nd": {"n": 60, "m": 144.24, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Pm": {"n": 61, "m": 145.0, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Sm": {"n": 62, "m": 150.36, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Eu": {"n": 63, "m": 151.96, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Gd": {"n": 64, "m": 157.25, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Tb": {"n": 65, "m": 158.93, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Dy": {"n": 66, "m": 162.5, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Ho": {"n": 67, "m": 164.93, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Er": {"n": 68, "m": 167.26, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Tm": {"n": 69, "m": 168.93, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Yb": {"n": 70, "m": 173.05, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Lu": {"n": 71, "m": 174.97, "cat": "Lantanídeo", "cor": "#ec4899"},
    "Hf": {"n": 72, "m": 178.49, "cat": "Transição", "cor": "#ef4444"},
    "Ta": {"n": 73, "m": 180.95, "cat": "Transição", "cor": "#ef4444"},
    "W": {"n": 74, "m": 183.84, "cat": "Transição", "cor": "#ef4444"},
    "Re": {"n": 75, "m": 186.21, "cat": "Transição", "cor": "#ef4444"},
    "Os": {"n": 76, "m": 190.23, "cat": "Transição", "cor": "#ef4444"},
    "Ir": {"n": 77, "m": 192.22, "cat": "Transição", "cor": "#ef4444"},
    "Pt": {"n": 78, "m": 195.08, "cat": "Transição", "cor": "#ef4444"},
    "Au": {"n": 79, "m": 196.97, "cat": "Transição", "cor": "#ef4444"},
    "Hg": {"n": 80, "m": 200.59, "cat": "Transição", "cor": "#ef4444"},
    "Tl": {"n": 81, "m": 204.38, "cat": "Outro metal", "cor": "#6b7280"},
    "Pb": {"n": 82, "m": 207.2, "cat": "Outro metal", "cor": "#6b7280"},
    "Bi": {"n": 83, "m": 208.98, "cat": "Outro metal", "cor": "#6b7280"},
    "Po": {"n": 84, "m": 209.0, "cat": "Semimetal", "cor": "#06b6d4"},
    "At": {"n": 85, "m": 210.0, "cat": "Halogênio", "cor": "#f43f5e"},
    "Rn": {"n": 86, "m": 222.0, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Fr": {"n": 87, "m": 223.0, "cat": "Alcalino", "cor": "#f59e0b"},
    "Ra": {"n": 88, "m": 226.0, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "Ac": {"n": 89, "m": 227.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Th": {"n": 90, "m": 232.04, "cat": "Actinídeo", "cor": "#f97316"},
    "Pa": {"n": 91, "m": 231.04, "cat": "Actinídeo", "cor": "#f97316"},
    "U": {"n": 92, "m": 238.03, "cat": "Actinídeo", "cor": "#f97316"},
    "Np": {"n": 93, "m": 237.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Pu": {"n": 94, "m": 244.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Am": {"n": 95, "m": 243.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Cm": {"n": 96, "m": 247.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Bk": {"n": 97, "m": 247.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Cf": {"n": 98, "m": 251.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Es": {"n": 99, "m": 252.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Fm": {"n": 100, "m": 257.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Md": {"n": 101, "m": 258.0, "cat": "Actinídeo", "cor": "#f97316"},
    "No": {"n": 102, "m": 259.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Lr": {"n": 103, "m": 266.0, "cat": "Actinídeo", "cor": "#f97316"},
    "Rf": {"n": 104, "m": 267.0, "cat": "Transição", "cor": "#ef4444"},
    "Db": {"n": 105, "m": 268.0, "cat": "Transição", "cor": "#ef4444"},
    "Sg": {"n": 106, "m": 271.0, "cat": "Transição", "cor": "#ef4444"},
    "Bh": {"n": 107, "m": 270.0, "cat": "Transição", "cor": "#ef4444"},
    "Hs": {"n": 108, "m": 277.0, "cat": "Transição", "cor": "#ef4444"},
    "Mt": {"n": 109, "m": 278.0, "cat": "Transição", "cor": "#ef4444"},
    "Ds": {"n": 110, "m": 281.0, "cat": "Transição", "cor": "#ef4444"},
    "Rg": {"n": 111, "m": 282.0, "cat": "Transição", "cor": "#ef4444"},
    "Cn": {"n": 112, "m": 285.0, "cat": "Transição", "cor": "#ef4444"},
    "Nh": {"n": 113, "m": 286.0, "cat": "Outro metal", "cor": "#6b7280"},
    "Fl": {"n": 114, "m": 289.0, "cat": "Outro metal", "cor": "#6b7280"},
    "Mc": {"n": 115, "m": 290.0, "cat": "Outro metal", "cor": "#6b7280"},
    "Lv": {"n": 116, "m": 293.0, "cat": "Outro metal", "cor": "#6b7280"},
    "Ts": {"n": 117, "m": 294.0, "cat": "Halogênio", "cor": "#f43f5e"},
    "Og": {"n": 118, "m": 294.0, "cat": "Gás Nobre", "cor": "#8b5cf6"},
    "Uue": {"n": 119, "m": 315.0, "cat": "Alcalino", "cor": "#f59e0b"},
    "Ubn": {"n": 120, "m": 320.0, "cat": "Alcalino-terroso", "cor": "#10b981"},
    "Ubu": {"n": 121, "m": 326.0, "cat": "Superpesado", "cor": "#7c3aed"},
    "Ubb": {"n": 122, "m": 328.0, "cat": "Superpesado", "cor": "#7c3aed"},
    "Ubt": {"n": 123, "m": 330.0, "cat": "Superpesado", "cor": "#7c3aed"},
    "Ubq": {"n": 124, "m": 332.0, "cat": "Superpesado", "cor": "#7c3aed"},
    "Ubp": {"n": 125, "m": 334.0, "cat": "Superpesado", "cor": "#7c3aed"},
    "Ubh": {"n": 126, "m": 336.0, "cat": "Superpesado", "cor": "#7c3aed"},
    "Ubs": {"n": 127, "m": 338.0, "cat": "Superpesado", "cor": "#7c3aed"}
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
        df_kps = pd.DataFrame([["AgCl", "1,6 x 10⁻¹⁰"], ["BaSO₄", "1,1 x 10⁻¹⁰"], ["CaCO₃", "3,36 x 10⁻⁹"],["PbBr₂", "7,9 x 10⁻⁵"],["CuBr", "4,2 x 10⁻⁸"],["AgBr", "7,7 x 10⁻¹³"],["Al(OH)₃", "1,1 x 10⁻³³"],["Fe(OH)₃", "4 x 10⁻³⁸"],["Mg(OH)₂", "1,8 x 10⁻¹¹"]], columns=["Fórmula", "Kps"])
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
    st.subheader("Comparação de Curvas de Solubilidade")
    
    if 'sais' not in st.session_state:
        st.session_state.sais = [
            {"nome": "KNO3", "temp": "0, 20, 40, 60, 80", "sol": "13, 32, 64, 110, 169", "cor": "#10b981"},
            {"nome": "NaCl", "temp": "0, 20, 40, 60, 80", "sol": "35, 36, 37, 38, 39", "cor": "#3b82f6"}
        ]

    c1, c2 = st.columns([1, 2])

    with c1:
        st.markdown("### 🛠️ Configurar Sais")
        df_sais = st.data_editor(st.session_state.sais, num_rows="dynamic", key="editor_sais")
        
        # Botão para baixar os dados em CSV
        import pandas as pd
        csv = pd.DataFrame(df_sais).to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Baixar Dados (CSV)",
            data=csv,
            file_name='curvas_solubilidade.csv',
            mime='text/csv',
        )

    with c2:
        try:
            fig = go.Figure()
            for sal in df_sais:
                x = np.array([float(i) for i in sal["temp"].split(",")])
                y = np.array([float(i) for i in sal["sol"].split(",")])
                x_smooth = np.linspace(x.min(), x.max(), 300)
                y_smooth = make_interp_spline(x, y, k=2)(x_smooth)
                
                fig.add_trace(go.Scatter(
                    x=x_smooth, y=y_smooth, name=sal["nome"],
                    line=dict(color=sal["cor"], width=3)
                ))

            fig.update_layout(
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                font=dict(color='white'),
                xaxis=dict(title="Temperatura (°C)", gridcolor='rgba(255,255,255,0.1)'),
                yaxis=dict(title="g/100g H2O", gridcolor='rgba(255,255,255,0.1)'),
                legend=dict(bgcolor='rgba(0,0,0,0)'),
                hovermode="x unified"
            )

            # Configuração para permitir o download da imagem pelo menu do gráfico
            config = {
                'toImageButtonOptions': {
                    'format': 'png', # ou 'jpeg', 'svg', 'pdf'
                    'filename': 'grafico_solubilidade',
                    'height': 500,
                    'width': 700,
                    'scale': 1 # Resolução da imagem
                }
            }

            st.plotly_chart(fig, use_container_width=True, config=config)
            st.caption("📸 Use a câmera no canto superior direito do gráfico para baixar como PNG.")

        except Exception as e:
            st.warning("Verifique o formato dos dados (ex: 10, 20, 30)")

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
