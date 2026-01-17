import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import plotly.express as px
import plotly.graph_objects as go

# -------------------------------
# Page Configuration
# -------------------------------
st.set_page_config(
    page_title="Gut Microbiome Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# -------------------------------
# Custom CSS for Better Styling
# -------------------------------
st.markdown("""
<style>
    /* Metric container */
    div[data-testid="metric-container"] {
        background-color: white;
        padding: 14px;
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.12);
        margin: 8px 6px;
    }
    div[data-testid="stMetricValue"] {
        font-size: 1.45rem;
        color: #2E86C1;
    }
    /* Page background */
    .main {
        background-color: #fbfbfb;
    }
    /* Sidebar styling */
    .css-1lcbmhc, [data-testid="stSidebar"] {
        background-color: #f0f2f6 !important;
    }
    
    /* Title and header text visibility */
    h1, h2, h3, h4, h5, h6 {
        color: #1a1a1a !important;
    }
    /* Body text */
    body, p, div, span {
        color: #333333 !important;
    }
    
    /* All text elements */
    * {
        color: #333333 !important;
    }
    
    /* Selectbox styling - target the select element itself */
    [data-baseweb="select"] {
        background-color: white !important;
    }
    [data-baseweb="select"] span, 
    [data-baseweb="select"] div {
        color: #1a1a1a !important;
    }
    
    /* Radio button text */
    .stRadio label, .stRadio span {
        color: #1a1a1a !important;
    }
    
    /* Checkbox text */
    .stCheckbox label, .stCheckbox span {
        color: #1a1a1a !important;
    }
    
    /* Slider label and text */
    .stSlider label, .stSlider span {
        color: #1a1a1a !important;
    }
    
    /* Selectbox label */
    .stSelectbox label {
        color: #1a1a1a !important;
    }
    
    /* Dropdown menu items text */
    [role="option"], 
    [role="listbox"] {
        color: #1a1a1a !important;
    }
    
    /* Navigation radio buttons */
    [data-testid="stRadio"] label {
        color: #1a1a1a !important;
    }
    
    /* All labels across the app */
    label {
        color: #1a1a1a !important;
    }
    
    /* DROPDOWN SELECTBOX - Most important styles */
    [data-baseweb="select"] > div {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    [data-baseweb="select"] input {
        color: #1a1a1a !important;
    }
    
    [data-baseweb="select"] span,
    [data-baseweb="select"] div,
    [data-baseweb="select"] p {
        color: #1a1a1a !important;
    }
    
    /* Popover menu items for selectbox - CRITICAL */
    [role="option"] {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    [role="option"] > div,
    [role="option"] span,
    [role="option"] p {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    /* Listbox styling with white background */
    [role="listbox"] {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    [role="listbox"] span,
    [role="listbox"] div,
    [role="listbox"] p {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    /* Popover for dropdown - make it white */
    [data-testid="stPopoverContent"],
    [class*="baseweb-popover"] {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    [data-testid="stPopoverContent"] span,
    [data-testid="stPopoverContent"] div,
    [data-testid="stPopoverContent"] p,
    [class*="baseweb-popover"] span,
    [class*="baseweb-popover"] div,
    [class*="baseweb-popover"] p {
        background-color: white !important;
        color: #1a1a1a !important;
    }
    
    /* Smaller text for footers */
    .small {
        font-size:12px;
        color:#555;
    }
    /* Improve table header visibility */
    .stDataFrame thead th {
        background-color: #f1f5f9;
        color: #1a1a1a !important;
    }
</style>
""", unsafe_allow_html=True)

# -------------------------------
# Title
# -------------------------------
st.title("🧬 Gut Microbiome Analysis Dashboard")
st.markdown("### Dynamic Nutrient Bioavailability Prediction Model")

# -------------------------------
# Load Data
# -------------------------------
@st.cache_data
def load_data(path="india_species_abundance_clean.tsv"):
    try:
        otu = pd.read_csv(path, sep="\t", index_col=0)
        # Ensure numeric and drop fully empty rows
        otu = otu.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        return otu
    except Exception as e:
        st.error("❌ Error loading OTU data. Make sure 'india_species_abundance_clean.tsv' is in the directory.")
        st.stop()

otu = load_data()
baseline_abundance = otu.mean(axis=1)

# -------------------------------
# Genus Mapping & Traits
# -------------------------------
def get_genus(species):
    # Extract genus from species name and convert to lowercase for matching
    if isinstance(species, str):
        genus = species.split("_")[0]
        return genus.lower()  # Convert to lowercase for consistent matching
    return "unknown"

species_list = baseline_abundance.index.tolist()
genus_map = {sp: get_genus(sp) for sp in species_list}

genus_traits = {
    "lactobacillus": {"SCFA":0.9,"pH_reduction":0.9,"Barrier_support":0.8,"Vitamin_Biosynthesis":0.6,"Siderophore":0.0},
    "bifidobacterium": {"SCFA":0.8,"pH_reduction":0.7,"Barrier_support":0.8,"Vitamin_Biosynthesis":0.7,"Siderophore":0.0},
    "faecalibacterium": {"SCFA":0.95,"pH_reduction":0.9,"Barrier_support":0.9,"Vitamin_Biosynthesis":0.4,"Siderophore":0.0},
    "roseburia": {"SCFA":0.9,"pH_reduction":0.8,"Barrier_support":0.7,"Vitamin_Biosynthesis":0.2,"Siderophore":0.0},
    "bacteroides": {"SCFA":0.6,"pH_reduction":0.3,"Barrier_support":0.5,"Vitamin_Biosynthesis":0.5,"Siderophore":0.1},
    "prevotella": {"SCFA":0.7,"pH_reduction":0.4,"Barrier_support":0.6,"Vitamin_Biosynthesis":0.3,"Siderophore":0.0},
    "streptococcus": {"SCFA":0.3,"pH_reduction":0.5,"Barrier_support":0.3,"Vitamin_Biosynthesis":0.3,"Siderophore":0.0},
    "clostridium": {"SCFA":0.8,"pH_reduction":0.6,"Barrier_support":0.4,"Vitamin_Biosynthesis":0.2,"Siderophore":0.0},
    "ruminococcus": {"SCFA":0.85,"pH_reduction":0.7,"Barrier_support":0.6,"Vitamin_Biosynthesis":0.3,"Siderophore":0.0},
    "escherichia": {"SCFA":0.1,"pH_reduction":0.0,"Barrier_support":0.1,"Vitamin_Biosynthesis":0.1,"Siderophore":0.9},
    "enterobacter": {"SCFA":0.1,"pH_reduction":0.0,"Barrier_support":0.1,"Vitamin_Biosynthesis":0.1,"Siderophore":0.8},
    "eubacterium": {"SCFA":0.5,"pH_reduction":0.4,"Barrier_support":0.3,"Vitamin_Biosynthesis":0.2,"Siderophore":0.0},
    "haemophilus": {"SCFA":0.2,"pH_reduction":0.2,"Barrier_support":0.2,"Vitamin_Biosynthesis":0.2,"Siderophore":0.3},
    "megasphaera": {"SCFA":0.7,"pH_reduction":0.5,"Barrier_support":0.4,"Vitamin_Biosynthesis":0.1,"Siderophore":0.0},
    "parasutterella": {"SCFA":0.3,"pH_reduction":0.2,"Barrier_support":0.2,"Vitamin_Biosynthesis":0.2,"Siderophore":0.2},
    "gemmiger": {"SCFA":0.6,"pH_reduction":0.4,"Barrier_support":0.5,"Vitamin_Biosynthesis":0.3,"Siderophore":0.0}
}

# Build species-level trait matrix (species x traits)
default_trait = {"SCFA":0.2,"pH_reduction":0.1,"Barrier_support":0.1,"Vitamin_Biosynthesis":0.1,"Siderophore":0.2}
all_traits = []
for sp in species_list:
    genus = genus_map.get(sp, None)
    row = genus_traits.get(genus, default_trait).copy()
    row["Species"] = sp
    all_traits.append(row)
traits = pd.DataFrame(all_traits).set_index("Species")

# -------------------------------
# Nutrient model (configurable)
# -------------------------------
nutrients = {
    "Iron": {"SCFA":0.4,"pH_reduction":0.3,"Barrier_support":0.2,"Siderophore":-0.6},
    "Vitamin_B12": {"Vitamin_Biosynthesis":0.7,"Barrier_support":0.2},
    "Folate": {"Vitamin_Biosynthesis":0.6},
    "Calcium": {"SCFA":0.5,"pH_reduction":0.4},
    "Magnesium": {"SCFA":0.4},
    "Zinc": {"Barrier_support":0.4,"Siderophore":-0.5}
}

def absorption_score(abundance, nutrient, normalize=True):
    """
    Compute absorption score for a nutrient given species abundance (pd.Series or dict).
    If normalize=True, divide by total abundance to get per-unit effect.
    """
    if isinstance(abundance, pd.Series):
        abundance = abundance.to_dict()
    score = 0.0
    for sp, ab in abundance.items():
        if sp in traits.index and ab is not None and ab > 0:
            for trait, w in nutrients[nutrient].items():
                if trait in traits.columns:
                    score += ab * traits.loc[sp, trait] * w
    if normalize:
        total_ab = sum([v for v in abundance.values() if v is not None])
        if total_ab > 0:
            return score / total_ab
        return 0.0
    return score

# -------------------------------
# Sidebar Navigation
# -------------------------------
st.sidebar.markdown("## 📊 Navigation")
page = st.sidebar.radio(
    "Go to:",
    ["Overview", "OTU Analysis", "Trait Coverage", "Nutrient Absorption", "Simulations", "About"],
    index=0
)

# -------------------------------
# Overview Page
# -------------------------------
if page == "Overview":
    st.header("🧬 Project Overview")

    with st.container():
        c1, c2, c3 = st.columns(3)
        c1.metric("🦠 Species Detected", len(baseline_abundance))
        c2.metric("📊 Samples Analyzed", otu.shape[1])
        c3.metric("🔬 Nutrient Models", len(nutrients))

    st.markdown("---")
    st.subheader("📋 About This Dashboard")
    st.write("""
    This microbiome analysis platform integrates:
    - **OTU Composition Analysis**
    - **Nutrient Bioavailability Modeling**
    - **Interactive Simulations**
    - **Functional Trait Mapping**
    """)

    st.markdown("---")
    st.subheader("🌟 Top 10 Most Abundant Species")
    top_species = baseline_abundance.sort_values(ascending=False).head(10)
    fig = px.bar(
        x=top_species.values,
        y=top_species.index,
        orientation='h',
        labels={"x": "Mean Relative Abundance", "y": "Species"},
        color=top_species.values,
        color_continuous_scale="Viridis"
    )
    fig.update_layout(height=520, margin=dict(l=40,r=40,t=60,b=40))
    st.plotly_chart(fig, use_container_width=True)

# -------------------------------
# OTU Analysis Page
# -------------------------------
elif page == "OTU Analysis":
    st.header("OTU Abundance Analysis")

    st.subheader("Abundance Distribution")
    col1, col2 = st.columns([1,1])
    with col1:
        fig = px.histogram(
            baseline_abundance,
            nbins=60,
            title="Distribution of Mean Species Abundance",
            labels={"value": "Mean Relative Abundance"},
            color_discrete_sequence=["#636EFA"]
        )
        fig.update_layout(height=420, margin=dict(l=40,r=40,t=60,b=40))
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Cumulative Abundance")
        sorted_abund = baseline_abundance.sort_values(ascending=False)
        cumsum = np.cumsum(sorted_abund.values)
        cumsum_pct = (cumsum / cumsum[-1]) * 100
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(
            y=cumsum_pct,
            mode='lines+markers',
            fill='tozeroy',
            name='Cumulative %',
            line=dict(color='#00CC96')
        ))
        fig2.update_layout(title="Cumulative Species Abundance", xaxis_title="Species (ranked)", yaxis_title="Cumulative Abundance (%)", height=420, margin=dict(l=40,r=40,t=60,b=40))
        st.plotly_chart(fig2, use_container_width=True)

    st.markdown("---")
    st.subheader("Species Data Table")
    display_df = pd.DataFrame({
        "Species": baseline_abundance.index,
        "Mean Abundance": baseline_abundance.values,
        "Genus": [genus_map.get(sp, "Unknown") for sp in baseline_abundance.index]
    }).sort_values("Mean Abundance", ascending=False).reset_index(drop=True)
    st.dataframe(display_df.style.background_gradient(subset=["Mean Abundance"], cmap="Blues"), use_container_width=True, height=420)

# -------------------------------
# Trait Coverage Page
# -------------------------------
elif page == "Trait Coverage":
    st.header("Trait Coverage & Diagnostics")
    st.write("Visualize how traits are distributed across species and check coverage.")

    # Heatmap of traits (transpose for readability)
    trait_matrix = traits.copy()
    if trait_matrix.shape[0] > 0:
        fig = px.imshow(
            trait_matrix.T,
            labels=dict(x="Species", y="Trait", color="Value"),
            aspect="auto",
            color_continuous_scale="Viridis"
        )
        fig.update_layout(height=520, margin=dict(l=40,r=40,t=60,b=40))
        st.plotly_chart(fig, use_container_width=True)

    st.markdown("---")
    st.subheader("Top Species by Trait (example: SCFA)")
    top_scfa = traits["SCFA"].sort_values(ascending=False).head(15)
    fig2 = px.bar(x=top_scfa.values, y=top_scfa.index, orientation='h', labels={"x":"SCFA trait value","y":"Species"}, color=top_scfa.values, color_continuous_scale="Blues")
    fig2.update_layout(height=420, margin=dict(l=40,r=40,t=60,b=40))
    st.plotly_chart(fig2, use_container_width=True)

# -------------------------------
# Nutrient Absorption Page
# -------------------------------
elif page == "Nutrient Absorption":
    st.header("Nutrient Bioavailability Analysis")

    # Baseline nutrient scores (normalized)
    nutrient_scores = {n: absorption_score(baseline_abundance, n, normalize=True) for n in nutrients.keys()}

    left, right = st.columns([1,2])
    with left:
        st.subheader("Baseline Absorption Scores")
        for nutrient, score in sorted(nutrient_scores.items(), key=lambda x: x[1], reverse=True):
            st.metric(nutrient, f"{score:.4f}")

    with right:
        st.subheader("Nutrient Comparison")
        fig = px.bar(
            x=list(nutrient_scores.keys()),
            y=list(nutrient_scores.values()),
            title="Baseline Nutrient Absorption Scores (normalized)",
            labels={"x":"Nutrient","y":"Absorption Score"},
            color=list(nutrient_scores.values()),
            color_continuous_scale="RdYlGn"
        )
        fig.update_layout(height=520, margin=dict(l=40,r=40,t=60,b=40))
        st.plotly_chart(fig, use_container_width=True)

    st.markdown("---")
    st.subheader("Nutrient-Trait Relationships")
    selected_nutrient = st.selectbox("Select Nutrient:", list(nutrients.keys()))

    # Compute per-species contribution (trait-weight only, independent of abundance)
    contrib = {}
    for sp in traits.index:
        score = 0.0
        for trait, weight in nutrients[selected_nutrient].items():
            if trait in traits.columns:
                score += traits.loc[sp, trait] * weight
        # Store all contributions, even if small
        contrib[sp] = max(score, 0)  # Use 0 if negative
    
    contrib_df = pd.DataFrame(list(contrib.items()), columns=["Species","Contribution"])
    contrib_df = contrib_df.sort_values("Contribution", ascending=False)
    # Show top species with contributions, filter only if exactly 0
    contrib_df = contrib_df[contrib_df["Contribution"] > 0].head(15)
    
    if len(contrib_df) > 0:
        fig3 = px.bar(contrib_df, x="Contribution", y="Species", orientation='h', 
                      color="Contribution", color_continuous_scale="RdBu")
        fig3.update_layout(height=520, margin=dict(l=40,r=40,t=60,b=40))
        st.plotly_chart(fig3, use_container_width=True)
    else:
        st.warning("No species with positive contributions for this nutrient.")

# -------------------------------
# Simulations Page
# -------------------------------
elif page == "Simulations":
    st.header("Microbiome Perturbation Simulator")
    st.write("Simulate the impact of adding or removing specific bacterial species on nutrient absorption (normalized scores).")

    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    with col1:
        selected_bacteria = st.selectbox("Select Bacteria:", sorted(list(set(baseline_abundance.index).intersection(traits.index))))
    with col2:
        selected_nutrient = st.selectbox("Select Nutrient:", list(nutrients.keys()), key="sim_nutrient")
    with col3:
        delta = st.slider("Δ Abundance (additive):", min_value=-0.2, max_value=0.5, value=0.05, step=0.01)

    def simulate_addition(abundance, microbe, delta):
        new = abundance.copy()
        if microbe in new:
            new[microbe] = max(0.0, new[microbe] + delta)
        else:
            new[microbe] = max(0.0, delta)
        return new

    baseline_score = absorption_score(baseline_abundance, selected_nutrient, normalize=True)
    new_abundance = simulate_addition(baseline_abundance.to_dict(), selected_bacteria, delta)
    new_score = absorption_score(new_abundance, selected_nutrient, normalize=True)
    change = new_score - baseline_score
    pct_change = (change / baseline_score * 100) if baseline_score != 0 else np.nan

    st.markdown("---")
    m1, m2, m3 = st.columns(3)
    m1.metric("Baseline Score", f"{baseline_score:.4f}")
    m2.metric("After Perturbation", f"{new_score:.4f}")
    
    # Color-coded impact metric: green for positive, red for negative
    impact_color = "🟢" if change >= 0 else "🔴"
    m3.metric("Change (Δ)", f"{change:+.4f}", delta=change, delta_color="inverse")
    
    # Display colored impact summary
    if change > 0:
        st.success(f"✅ **Positive Impact**: +{change:.4f} ({pct_change:+.1f}%)")
    elif change < 0:
        st.error(f"❌ **Negative Impact**: {change:.4f} ({pct_change:+.1f}%)")
    else:
        st.info(f"ℹ️ **No Change**: {change:.4f} (0.0%)")

    st.markdown("---")
    chart_col, summary_col = st.columns([1.6, 1])
    with chart_col:
        comparison_data = pd.DataFrame({
            "Condition": ["Baseline", "Perturbed"],
            "Score": [baseline_score, new_score]
        })
        fig = px.bar(comparison_data, x="Condition", y="Score", color="Score", color_continuous_scale="RdYlGn", text="Score")
        fig.update_traces(texttemplate="%{text:.4f}", textposition='outside')
        fig.update_layout(height=480, margin=dict(l=40,r=40,t=60,b=40), showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    with summary_col:
        st.info(f"""
**Simulation Summary**

**Bacteria:** {selected_bacteria}

**Δ Abundance:** {delta:+.2f}

**Nutrient:** {selected_nutrient}

**Impact:** {change:+.4f} ({pct_change:+.1f}%)
""")

# -------------------------------
# About Page
# -------------------------------
elif page == "About":
    st.header("About This Project")
    st.subheader("Overview")
    st.write("""
    This microbiome analysis dashboard models the relationship between gut bacterial composition 
    and nutrient bioavailability using functional traits and interaction networks.
    """)
    st.subheader("Methodology")
    st.write("""
    **Data Processing:**
    - OTU tables normalized to relative abundance
    - Spearman correlation networks for co-occurrence patterns

    **Functional Traits:**
    - Genus-level priors based on published literature
    - SCFA production, pH reduction, barrier support, vitamin biosynthesis, siderophore production

    **Nutrient Model:**
    - Trait-to-nutrient mapping through trait weight relationships
    - Absorption scores reflect potential nutrient bioavailability

    **Simulations:**
    - Real-time response to microbiome perturbations
    """)
    st.subheader("Data Requirements")
    st.write("- `india_species_abundance_clean.tsv` - OTU abundance table (species × samples)")
    st.subheader("Built With")
    st.write("🐍 Python | 📊 Streamlit | 📈 Plotly | 🔬 SciPy")

# -------------------------------
# Footer
# -------------------------------
st.markdown("---")
st.markdown('<div class="small">🧬 Gut Microbiome Analysis Dashboard | Built with Streamlit</div>', unsafe_allow_html=True)
