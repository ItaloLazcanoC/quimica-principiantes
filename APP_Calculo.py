import math
import re
import streamlit as st

try:
    import periodictable as pt
except ImportError:
    pt = None

AVOGADRO = 6.02214076e23  # mol^-1

# ===========================
# Utilidades
# ===========================

class FormulaParseError(Exception):
    pass

def get_atomic_weight(symbol: str) -> float:
    if pt is None:
        raise RuntimeError("Falta la librería 'periodictable'. Instálala con: pip install periodictable")
    if symbol == "D":  # deuterio (opcional)
        return 2.01410177812
    try:
        el = getattr(pt, symbol)
        return float(el.mass)
    except Exception:
        raise FormulaParseError(f"Símbolo químico desconocido: '{symbol}'. Verifica mayúsculas/minúsculas (ej: 'Na', no 'NA').")

_token_pat = re.compile(r"""
    (?P<dot>[\.\·])|
    (?P<lparen>[\(\[\{])|
    (?P<rparen>[\)\]\}])|
    (?P<number>\d+)|
    (?P<elem>[A-Z][a-z]?)|
    (?P<space>\s+)
""", re.VERBOSE)

def tokenize(formula: str):
    pos = 0
    while pos < len(formula):
        m = _token_pat.match(formula, pos)
        if not m:
            raise FormulaParseError(f"Revisa la fórmula: hay un símbolo no reconocido cerca de '{formula[pos:pos+5]}'.")
        pos = m.end()
        kind = m.lastgroup
        val = m.group(m.lastgroup)
        if kind != "space":
            yield (kind, val)

def multiply_counts(counts, k):
    return {el: cnt * k for el, cnt in counts.items()}

def merge_counts(a, b):
    out = dict(a)
    for el, cnt in b.items():
        out[el] = out.get(el, 0) + cnt
    return out

def parse_formula(formula: str) -> dict:
    tokens = list(tokenize(formula))
    i = 0
    def parse_unit():
        nonlocal i
        counts = {}
        while i < len(tokens):
            kind, val = tokens[i]
            if kind == "rparen":
                break
            if kind == "dot":
                i += 1
                right = parse_unit()
                counts.update(merge_counts(counts, right))
                break
            elif kind == "elem":
                i += 1
                el = val
                mult = 1
                if i < len(tokens) and tokens[i][0] == "number":
                    mult = int(tokens[i][1]); i += 1
                counts[el] = counts.get(el, 0) + mult
            elif kind == "lparen":
                i += 1
                inner = parse_unit()
                if i >= len(tokens) or tokens[i][0] != "rparen":
                    raise FormulaParseError("Falta cerrar un paréntesis.")
                i += 1
                mult = 1
                if i < len(tokens) and tokens[i][0] == "number":
                    mult = int(tokens[i][1]); i += 1
                inner = multiply_counts(inner, mult)
                counts = merge_counts(counts, inner)
            elif kind == "number":
                mult = int(val); i += 1
                if i >= len(tokens) or tokens[i][0] not in ("elem", "lparen"):
                    raise FormulaParseError("Hay un número sin un elemento o grupo después.")
                if tokens[i][0] == "elem":
                    _, el = tokens[i]; i += 1
                    cnt = 1
                    if i < len(tokens) and tokens[i][0] == "number":
                        cnt = int(tokens[i][1]); i += 1
                    counts[el] = counts.get(el, 0) + mult * cnt
                else:
                    i += 1
                    inner = parse_unit()
                    if i >= len(tokens) or tokens[i][0] != "rparen":
                        raise FormulaParseError("Falta cerrar un paréntesis tras el número.")
                    i += 1
                    cnt = 1
                    if i < len(tokens) and tokens[i][0] == "number":
                        cnt = int(tokens[i][1]); i += 1
                    inner = multiply_counts(inner, mult * cnt)
                    counts = merge_counts(counts, inner)
            else:
                i += 1
        return counts
    total = parse_unit()
    if i < len(tokens) and tokens[i][0] == "rparen":
        raise FormulaParseError("Sobra un paréntesis de cierre.")
    return total

def molar_mass(formula: str) -> float:
    counts = parse_formula(formula)
    return sum(get_atomic_weight(el) * n for el, n in counts.items())

def percent_composition(formula: str):
    counts = parse_formula(formula)
    total_mm = molar_mass(formula)
    comp = {el: 100.0 * get_atomic_weight(el) * n / total_mm for el, n in counts.items()}
    return comp, total_mm

def empirical_formula_from_pairs(pairs, tolerance=0.05):
    total = sum(m for _, m in pairs)
    base100 = 95.0 <= total <= 105.0
    moles = []
    for sym, val in pairs:
        mm = get_atomic_weight(sym)
        grams = val if not base100 else val  # base 100 g si son %
        moles.append((sym, grams / mm))
    minmol = min(v for _, v in moles if v > 0)
    ratios = [(sym, v / minmol) for sym, v in moles]
    for k in [1, 2, 3, 4, 5, 6, 8, 10]:
        ok = True; rounded = []
        for sym, r in ratios:
            x = r * k; n = round(x)
            if n == 0 or abs(x - n) > tolerance:
                ok = False; break
            rounded.append((sym, n))
        if ok:
            counts = {}
            for sym, n in rounded:
                counts[sym] = counts.get(sym, 0) + n
            return format_formula(counts)
    counts = {sym: max(1, int(round(r))) for sym, r in ratios}
    return format_formula(counts)

def molecular_formula(empirical_formula: str, target_molar_mass: float):
    emp_counts = parse_formula(empirical_formula)
    emp_mm = sum(get_atomic_weight(el) * n for el, n in emp_counts.items())
    ratio = target_molar_mass / emp_mm
    k = int(round(ratio))
    if k <= 0 or abs(ratio - k) > 0.03:
        raise ValueError("La masa molar objetivo no es un múltiplo entero de la fórmula empírica (tolerancia ~3%).")
    for el in emp_counts:
        emp_counts[el] *= k
    return format_formula(emp_counts)

def format_formula(counts: dict) -> str:
    def order_key(item):
        el, _ = item
        if el == "C": return (0, el)
        if el == "H": return (1, el)
        return (2, el)
    parts = []
    for el, n in sorted(counts.items(), key=order_key):
        parts.append(f"{el}{'' if n==1 else n}")
    return "".join(parts)

# ===========================
# Contenidos Didácticos
# ===========================

GLOSARIO = {
    "Mol": "Es una 'cantidad' de materia. 1 mol = 6,022×10²³ partículas (átomos o moléculas).",
    "Número de Avogadro": "Un mol de una sustancia es igual a 6.022 × 10²³ unidades de esa sustancia (tal como átomos, moléculas, o iones). El número 6.022 × 10²³ se conoce como número de Avogadro o constante de Avogadro.",
    "Masa molar": "La masa molar (M) es la cantidad de masa que una sustancia contiene en un mol. Un mol se define como 6.022 x10 23 partículas",
    "Composición porcentual": "Qué porcentaje de la masa total aporta cada elemento del compuesto.",
    "Fórmula empírica": "Es una expresión que representa los átomos que forman un compuesto químico sin atender a su estructura. Es por tanto la representación más sencilla de un compuesto.",
    "Fórmula molecular": "Una fórmula molecular es una representación concisa de la composición de un compuesto químico, que especifica los tipos y números de átomos presentes en una molécula. Por ejemplo, la fórmula molecular del metanol es CH₄O, que indica un átomo de carbono, cuatro átomos de hidrógeno y uno de oxígeno."
}

EXAMPLES = {
    "Agua": "H2O",
    "Glucosa": "C6H12O6",
    "Sulfato de cobre pentahidratado": "CuSO4·5H2O",
    "Sulfato de amonio": "(NH4)2SO4",
    "Fosfato de calcio": "Ca3(PO4)2",
}

# ===========================
# UI
# ===========================

st.set_page_config(page_title="Química para principiantes", page_icon="🧪", layout="centered")

st.title("🧪 Química para principiantes: moles, masas y fórmulas")
st.markdown(
    "Esta herramienta **te guía paso a paso** para resolver ejercicios de:\n"
    "- Moles, Número de Avogadro, átomos/moléculas\n"
    "- Masa molar y **composición porcentual**\n"
    "- **Fórmula empírica** y **fórmula molecular**\n\n"
   
)

with st.sidebar:
    st.header("📚 Glosario rápido")
    for k, v in GLOSARIO.items():
        with st.expander(k):
            st.write(v)
    st.divider()
    st.subheader("🧪 Ejemplos de fórmulas")
    for name, f in EXAMPLES.items():
        if st.button(f"Usar {name}", key=f"ex_{name}"):
            st.session_state["last_formula"] = f
    st.caption("Consejo: las letras respetan mayúsculas/minúsculas (Na ≠ NA).")

tabs = st.tabs(["👣 Modo guiado", "🧮 Cálculos rápidos", "⚖️ Masa molar y %", "🧩 Empírica y molecular", "📘 Ejemplos resueltos"])

# ---------- Tab 1: Modo guiado ----------
with tabs[0]:
    st.subheader("👣 Paso a paso")
    step = st.radio("¿Qué necesitas?", [
        "¿Qué es un mol? (moles ↔ partículas)",
        "Gramos ↔ Moles (con una fórmula)",
        "Calcular masa molar y % composición",
        "Fórmula empírica desde % o masas",
        "Fórmula molecular desde empírica + masa molar"
    ])

    if step == "¿Qué es un mol? (moles ↔ partículas)":
        st.info("**Idea clave:** 1 mol = 6,022×10²³ partículas (Número de Avogadro).")
        op = st.selectbox("Elige conversión:", ["Moles → Partículas", "Partículas → Moles"])
        if op == "Moles → Partículas":
            n = st.number_input("¿Cuántos moles?", min_value=0.0, value=1.0)
            if st.button("Calcular", key="g1"):
                p = n * AVOGADRO
                st.latex(r"N_{\text{partículas}} = n \times N_A")
                st.success(f"Resultado: **{p:.3e} partículas**")
                st.caption("Multiplicamos los moles por el Número de Avogadro.")
        else:
            p = st.number_input("¿Cuántas partículas?", min_value=0.0, value=AVOGADRO)
            if st.button("Calcular", key="g2"):
                n = p / AVOGADRO
                st.latex(r"n = \dfrac{N_{\text{partículas}}}{N_A}")
                st.success(f"Resultado: **{n:.6f} mol**")
                st.caption("Dividimos las partículas por el Número de Avogadro.")

    if step == "Gramos ↔ Moles (con una fórmula)":
        formula = st.text_input("Escribe la fórmula química (ej: H2O, C6H12O6)", st.session_state.get("last_formula","H2O"))
        op2 = st.selectbox("Elige conversión:", ["Gramos → Moles", "Moles → Gramos"])
        if op2 == "Gramos → Moles":
            g = st.number_input("Masa (g)", min_value=0.0, value=18.0)
            if st.button("Calcular", key="g3"):
                try:
                    mm = molar_mass(formula)
                    n = g / mm
                    st.latex(r"n = \dfrac{m}{M}")
                    st.success(f"Masa molar de {formula}: **{mm:.5f} g/mol**\n\n"
                               f"{g:.4f} g → **{n:.6f} mol**")
                    st.caption("Dividimos la masa en gramos por la masa molar (g/mol).")
                except Exception as e:
                    st.error(str(e))
        else:
            n = st.number_input("Cantidad (mol)", min_value=0.0, value=1.0)
            if st.button("Calcular", key="g4"):
                try:
                    mm = molar_mass(formula)
                    g = n * mm
                    st.latex(r"m = n \times M")
                    st.success(f"Masa molar de {formula}: **{mm:.5f} g/mol**\n\n"
                               f"{n:.6f} mol → **{g:.5f} g**")
                    st.caption("Multiplicamos los moles por la masa molar.")
                except Exception as e:
                    st.error(str(e))

    if step == "Calcular masa molar y % composición":
        formula = st.text_input("Fórmula química (ej: CuSO4·5H2O)", st.session_state.get("last_formula","C6H12O6"), key="mmf")
        if st.button("Calcular", key="g5"):
            try:
                comp, mm = percent_composition(formula)
                st.latex(r"M = \sum (\text{masa atómica} \times \text{cantidad})")
                st.success(f"Masa molar de {formula}: **{mm:.5f} g/mol**")
                st.markdown("**Composición porcentual (aporte de cada elemento):**")
                for el, pct in sorted(comp.items()):
                    st.write(f"- {el}: **{pct:.4f}%**")
                st.caption("Cada porcentaje indica cuánta masa del compuesto corresponde a ese elemento.")
            except Exception as e:
                st.error(str(e))

    if step == "Fórmula empírica desde % o masas":
        st.info("Tip: si ingresas porcentajes que suman ~100, se asume **base 100 g**.")
        nrows = st.number_input("¿Cuántos elementos hay?", min_value=2, max_value=8, value=3)
        pairs = []
        for i in range(int(nrows)):
            c1, c2 = st.columns([1,2])
            sym = c1.text_input(f"Símbolo {i+1}", value="CHO"[i] if i<3 else "")
            val = c2.number_input(f"Valor {i+1} (masa en g o %)", min_value=0.0, value=[40.0,6.7,53.3][i] if i<3 else 0.0)
            if sym:
                pairs.append((sym.strip(), val))
        if st.button("Obtener fórmula empírica", key="g6"):
            try:
                emp = empirical_formula_from_pairs(pairs)
                st.latex(r"\text{Cálculo: convertir a moles, dividir por el menor y redondear a enteros simples.}")
                st.success(f"Fórmula empírica: **{emp}**")
                st.caption("Se convierten las masas/% a moles y se buscan las proporciones más simples.")
            except Exception as e:
                st.error(str(e))

    if step == "Fórmula molecular desde empírica + masa molar":
        emp_in = st.text_input("Fórmula empírica (ej: CH2O)", "CH2O")
        target_mm = st.number_input("Masa molar objetivo (g/mol)", min_value=0.0, value=180.16)
        if st.button("Calcular fórmula molecular", key="g7"):
            try:
                mf = molecular_formula(emp_in, target_mm)
                st.latex(r"\text{La fórmula molecular = } k \times \text{(fórmula empírica)}")
                st.success(f"Fórmula molecular: **{mf}**")
                st.caption("Buscamos un múltiplo entero k tal que k·M(empírica) ≈ Masa molar objetivo.")
            except Exception as e:
                st.error(str(e))

# ---------- Tab 2: Cálculos rápidos ----------
with tabs[1]:
    st.subheader("🧮 Conversor rápido")
    conv = st.selectbox("¿Qué quieres convertir?", [
        "Gramos → Moles",
        "Moles → Gramos",
        "Moles → Partículas",
        "Partículas → Moles",
    ])
    if "Gramos" in conv or "Moles" in conv:
        formula = st.text_input("Fórmula química", st.session_state.get("last_formula","H2O"), key="fquick")
    if conv == "Gramos → Moles":
        g = st.number_input("Masa (g)", min_value=0.0, value=18.0)
        if st.button("Calcular", key="q1"):
            try:
                mm = molar_mass(formula); n = g/mm
                st.success(f"M({formula})={mm:.5f} g/mol → **{n:.6f} mol**")
            except Exception as e:
                st.error(str(e))
    elif conv == "Moles → Gramos":
        n = st.number_input("Cantidad (mol)", min_value=0.0, value=1.0)
        if st.button("Calcular", key="q2"):
            try:
                mm = molar_mass(formula); g = n*mm
                st.success(f"M({formula})={mm:.5f} g/mol → **{g:.5f} g**")
            except Exception as e:
                st.error(str(e))
    elif conv == "Moles → Partículas":
        n = st.number_input("Cantidad (mol)", min_value=0.0, value=1.0)
        if st.button("Calcular", key="q3"):
            st.success(f"**{n*AVOGADRO:.3e} partículas**")
    elif conv == "Partículas → Moles":
        p = st.number_input("Número de partículas", min_value=0.0, value=AVOGADRO)
        if st.button("Calcular", key="q4"):
            st.success(f"**{p/AVOGADRO:.6f} mol**")

# ---------- Tab 3: Masa molar y % ----------
with tabs[2]:
    st.subheader("⚖️ Masa molar y composición")
    formula2 = st.text_input("Fórmula química (ej: C6H12O6, CuSO4·5H2O)", st.session_state.get("last_formula","C6H12O6"), key="f2")
    if st.button("Calcular", key="m1"):
        try:
            comp, mm = percent_composition(formula2)
            st.success(f"Masa molar de {formula2}: **{mm:.5f} g/mol**")
            st.markdown("**Composición porcentual:**")
            for el, pct in sorted(comp.items()):
                st.write(f"- {el}: **{pct:.4f}%**")
        except Exception as e:
            st.error(str(e))

# ---------- Tab 4: Empírica y molecular ----------
with tabs[3]:
    st.subheader("🧩 Fórmulas empírica y molecular")
    st.markdown("**Empírica desde % o masas**")
    nrows = st.number_input("Cantidad de elementos", min_value=2, max_value=8, value=3, key="nrows2")
    pairs = []
    for i in range(int(nrows)):
        c1, c2 = st.columns([1,2])
        sym = c1.text_input(f"Símbolo {i+1}", value="CHO"[i] if i<3 else "", key=f"s2_{i}")
        val = c2.number_input(f"Valor {i+1} (g o %)", min_value=0.0, value=[40.0,6.7,53.3][i] if i<3 else 0.0, key=f"v2_{i}")
        if sym:
            pairs.append((sym.strip(), val))
    if st.button("Calcular empírica", key="e1"):
        try:
            emp = empirical_formula_from_pairs(pairs)
            st.success(f"Empírica: **{emp}**")
        except Exception as e:
            st.error(str(e))

    st.markdown("---")
    st.markdown("**Molecular desde empírica + masa molar**")
    emp_in = st.text_input("Fórmula empírica", "CH2O", key="emp_in2")
    target_mm = st.number_input("Masa molar objetivo (g/mol)", min_value=0.0, value=180.156, key="tmm2")
    if st.button("Calcular molecular", key="e2"):
        try:
            mf = molecular_formula(emp_in, target_mm)
            st.success(f"Molecular: **{mf}**")
        except Exception as e:
            st.error(str(e))

# ---------- Tab 5: Ejemplos resueltos ----------
with tabs[4]:
    st.subheader("📘 Ejemplos paso a paso")
    with st.expander("1) ¿Cuántos moles hay en 36,0 g de agua (H2O)?"):
        try:
            mm = molar_mass("H2O")
            n = 36.0 / mm
            st.write(f"M(H2O) = {mm:.4f} g/mol → **{n:.4f} mol**")
            st.caption("Fórmula usada: n = m / M")
        except Exception as e:
            st.warning(str(e))

    with st.expander("2) ¿Cuántas moléculas hay en 0,250 mol de CO2?"):
        st.write(f"Moléculas = 0,250 × {AVOGADRO:.6e} = **{0.250*AVOGADRO:.3e}**")

    with st.expander("3) Masa molar y % en CuSO4·5H2O"):
        try:
            mm = molar_mass("CuSO4·5H2O")
            comp, _ = percent_composition("CuSO4·5H2O")
            st.write(f"M = **{mm:.4f} g/mol**")
            for el, pct in sorted(comp.items()):
                st.write(f"- {el}: {pct:.4f}%")
        except Exception as e:
            st.warning(str(e))

    with st.expander("4) Fórmula empírica desde %: C=40,00; H=6,71; O=53,29"):
        try:
            emp = empirical_formula_from_pairs([("C",40.00),("H",6.71),("O",53.29)])
            st.write(f"Empírica = **{emp}** (esperada: CH2O)")
        except Exception as e:
            st.warning(str(e))

    with st.expander("5) Fórmula molecular: empírica CH2O y masa molar 180,16 g/mol"):
        try:
            mf = molecular_formula("CH2O", 180.16)
            st.write(f"Molecular = **{mf}** (glucosa)")
        except Exception as e:
            st.warning(str(e))

st.divider()
st.caption("Diseñada para uso educativo inicial. Las masas atómicas provienen de la librería 'periodictable'.")
