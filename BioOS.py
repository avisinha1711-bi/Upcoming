"""
BIOLAB-OS: Biological Laboratory Operating System
A comprehensive OS for managing virology, genetics, and molecular biology labs
"""

import json
import sqlite3
import hashlib
import datetime
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
import threading
import queue
import time
import uuid
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

# ==================== DATABASE MODULE ====================

class BioDatabase:
    """Centralized database for all biological data"""
    
    def __init__(self, db_path="biolab.db"):
        self.conn = sqlite3.connect(db_path, check_same_thread=False)
        self.create_tables()
    
    def create_tables(self):
        """Initialize all database tables"""
        cursor = self.conn.cursor()
        
        # Users and authentication
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS users (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                username TEXT UNIQUE NOT NULL,
                password_hash TEXT NOT NULL,
                role TEXT NOT NULL,
                department TEXT,
                email TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Samples
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS samples (
                id TEXT PRIMARY KEY,
                sample_type TEXT NOT NULL,
                source TEXT,
                collection_date DATE,
                location TEXT,
                storage_conditions TEXT,
                volume_ul REAL,
                concentration_ng_ul REAL,
                purity_260_280 REAL,
                status TEXT DEFAULT 'active',
                created_by INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (created_by) REFERENCES users (id)
            )
        ''')
        
        # Viruses
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS viruses (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                family TEXT,
                genus TEXT,
                host_range TEXT,
                genome_type TEXT,
                genome_length INTEGER,
                sequence_file TEXT,
                virulence_factors TEXT,
                transmission TEXT,
                isolation_date DATE,
                lab_location TEXT,
                risk_level TEXT,
                created_by INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (created_by) REFERENCES users (id)
            )
        ''')
        
        # DNA Sequences
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sequences (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                sequence TEXT NOT NULL,
                length INTEGER,
                gc_content REAL,
                type TEXT,
                gene_name TEXT,
                organism TEXT,
                annotations TEXT,
                created_by INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (created_by) REFERENCES users (id)
            )
        ''')
        
        # PCR Experiments
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS pcr_experiments (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                target_gene TEXT,
                forward_primer TEXT,
                reverse_primer TEXT,
                annealing_temp REAL,
                cycle_count INTEGER,
                template_dna TEXT,
                result_file TEXT,
                ct_value REAL,
                amplification_efficiency REAL,
                status TEXT,
                created_by INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (created_by) REFERENCES users (id)
            )
        ''')
        
        # CRISPR Experiments
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS crispr_experiments (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                target_sequence TEXT,
                guide_rna TEXT,
                cas_type TEXT,
                cell_line TEXT,
                efficiency REAL,
                off_target_effects TEXT,
                validation_method TEXT,
                results TEXT,
                created_by INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (created_by) REFERENCES users (id)
            )
        ''')
        
        # Inventory
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS inventory (
                id TEXT PRIMARY KEY,
                item_name TEXT NOT NULL,
                catalog_number TEXT,
                supplier TEXT,
                quantity INTEGER,
                unit TEXT,
                storage_temp TEXT,
                expiration_date DATE,
                location TEXT,
                reorder_level INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Experiments
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS experiments (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                type TEXT,
                objective TEXT,
                protocol TEXT,
                materials_used TEXT,
                results TEXT,
                conclusion TEXT,
                status TEXT DEFAULT 'in_progress',
                created_by INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                completed_at TIMESTAMP,
                FOREIGN KEY (created_by) REFERENCES users (id)
            )
        ''')
        
        self.conn.commit()

# ==================== CORE MODULES ====================

class BioSecurityLevel(Enum):
    """Biosafety levels"""
    BSL1 = "BSL-1"
    BSL2 = "BSL-2"
    BSL3 = "BSL-3"
    BSL4 = "BSL-4"

@dataclass
class BioSample:
    """Biological sample container"""
    id: str
    sample_type: str  # DNA, RNA, Protein, Virus, Bacteria, Tissue
    source: str
    collection_date: datetime.date
    volume_ul: float
    concentration_ng_ul: Optional[float]
    purity_260_280: Optional[float]
    storage_temp: str
    location: str
    metadata: Dict

class VirusAnalyzer:
    """Virus analysis and characterization module"""
    
    @staticmethod
    def calculate_genome_characteristics(sequence: str) -> Dict:
        """Analyze viral genome"""
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
        at_content = 100 - gc_content
        
        return {
            'length': len(sequence),
            'gc_content': gc_content,
            'at_content': at_content,
            'gc_skew': (sequence.count('G') - sequence.count('C')) / (sequence.count('G') + sequence.count('C')) if (sequence.count('G') + sequence.count('C')) > 0 else 0,
            'codon_usage': VirusAnalyzer._calculate_codon_usage(sequence)
        }
    
    @staticmethod
    def predict_open_reading_frames(sequence: str, min_length=300):
        """Find potential ORFs in viral genome"""
        orfs = []
        for frame in range(3):
            for start in range(frame, len(sequence)-2, 3):
                if sequence[start:start+3] == "ATG":  # Start codon
                    for stop in range(start+3, len(sequence)-2, 3):
                        codon = sequence[stop:stop+3]
                        if codon in ["TAA", "TAG", "TGA"]:  # Stop codons
                            orf_length = stop - start
                            if orf_length >= min_length:
                                orfs.append({
                                    'start': start,
                                    'stop': stop,
                                    'length': orf_length,
                                    'frame': frame,
                                    'sequence': sequence[start:stop+3]
                                })
                            break
        return orfs
    
    @staticmethod
    def _calculate_codon_usage(sequence: str) -> Dict:
        """Calculate codon usage frequency"""
        codon_table = {}
        for i in range(0, len(sequence)-2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                codon_table[codon] = codon_table.get(codon, 0) + 1
        return codon_table

class DNASequenceAnalyzer:
    """DNA sequence analysis module"""
    
    @staticmethod
    def align_sequences(seq1: str, seq2: str) -> Tuple[str, str, float]:
        """Simple Needleman-Wunsch alignment"""
        # Simplified implementation
        gap_penalty = -1
        match_score = 1
        mismatch_score = -1
        
        # Create scoring matrix
        n, m = len(seq1), len(seq2)
        score = np.zeros((n+1, m+1))
        
        # Initialize
        for i in range(n+1):
            score[i][0] = i * gap_penalty
        for j in range(m+1):
            score[0][j] = j * gap_penalty
        
        # Fill matrix
        for i in range(1, n+1):
            for j in range(1, m+1):
                match = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
                delete = score[i-1][j] + gap_penalty
                insert = score[i][j-1] + gap_penalty
                score[i][j] = max(match, delete, insert)
        
        # Traceback
        align1, align2 = "", ""
        i, j = n, m
        
        while i > 0 or j > 0:
            if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif i > 0 and score[i][j] == score[i-1][j] + gap_penalty:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1
        
        identity = sum(1 for a, b in zip(align1, align2) if a == b) / len(align1) * 100
        return align1, align2, identity
    
    @staticmethod
    def find_restriction_sites(sequence: str) -> Dict:
        """Find restriction enzyme sites"""
        restriction_enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'NotI': 'GCGGCCGC',
            'XbaI': 'TCTAGA',
            'SpeI': 'ACTAGT'
        }
        
        sites = {}
        for enzyme, site in restriction_enzymes.items():
            positions = []
            for i in range(len(sequence) - len(site) + 1):
                if sequence[i:i+len(site)] == site:
                    positions.append(i)
            if positions:
                sites[enzyme] = positions
        
        return sites

class PCRSimulator:
    """PCR experiment simulation and optimization"""
    
    @staticmethod
    def simulate_pcr(sequence: str, forward_primer: str, reverse_primer: str, 
                    cycles: int = 30) -> Dict:
        """Simulate PCR amplification"""
        
        # Check primer binding
        forward_binding = sequence.find(forward_primer)
        reverse_binding = sequence.find(DNASequenceAnalyzer.reverse_complement(reverse_primer))
        
        if forward_binding == -1 or reverse_binding == -1:
            return {'success': False, 'error': 'Primers do not bind'}
        
        # Calculate product size
        product_size = reverse_binding - forward_binding + len(reverse_primer)
        
        # Simulate amplification
        initial_template = 1
        final_product = initial_template * (2 ** cycles)
        
        # Calculate efficiency
        efficiency = 0.95  # Typical efficiency
        
        return {
            'success': True,
            'product_size': product_size,
            'forward_binding_position': forward_binding,
            'reverse_binding_position': reverse_binding,
            'cycles': cycles,
            'estimated_product_yield': final_product,
            'efficiency': efficiency
        }
    
    @staticmethod
    def optimize_temperature(tm: float) -> Dict:
        """Optimize PCR temperatures"""
        return {
            'denaturation': 94.0,
            'annealing': tm - 5,  # 5°C below Tm
            'extension': 72.0,
            'gradient_start': tm - 10,
            'gradient_end': tm + 5
        }

class CRISPRDesigner:
    """CRISPR guide RNA design module"""
    
    @staticmethod
    def design_guide_rna(target_sequence: str, pam_sequence: str = "NGG") -> List[Dict]:
        """Design CRISPR guide RNAs for a target sequence"""
        guides = []
        
        for i in range(len(target_sequence) - 19):  # 20nt guide + 3nt PAM
            potential_guide = target_sequence[i:i+20]
            
            # Check for PAM site
            if pam_sequence == "NGG":
                pam = target_sequence[i+20:i+23]
                if pam[1:] == "GG":  # NGG PAM
                    guides.append({
                        'position': i,
                        'guide_sequence': potential_guide,
                        'pam_sequence': pam,
                        'gc_content': (potential_guide.count('G') + potential_guide.count('C')) / 20 * 100,
                        'specificity_score': CRISPRDesigner._calculate_specificity(potential_guide),
                        'off_target_risk': CRISPRDesigner._assess_off_target_risk(potential_guide)
                    })
        
        # Sort by specificity score
        guides.sort(key=lambda x: x['specificity_score'], reverse=True)
        return guides[:10]  # Return top 10
    
    @staticmethod
    def _calculate_specificity(guide_sequence: str) -> float:
        """Calculate guide RNA specificity score"""
        gc_content = (guide_sequence.count('G') + guide_sequence.count('C')) / len(guide_sequence)
        
        # Penalize extreme GC content
        gc_score = 1 - abs(0.5 - gc_content)
        
        # Score based on sequence complexity
        complexity = len(set(guide_sequence)) / len(guide_sequence)
        
        return (gc_score * 0.6 + complexity * 0.4) * 100
    
    @staticmethod
    def _assess_off_target_risk(guide_sequence: str) -> str:
        """Assess off-target risk level"""
        # Simple heuristic based on GC content and repeats
        gc_content = (guide_sequence.count('G') + guide_sequence.count('C')) / len(guide_sequence)
        
        if gc_content < 0.3 or gc_content > 0.7:
            return "High"
        elif 'AAAA' in guide_sequence or 'CCCC' in guide_sequence or 'GGGG' in guide_sequence or 'TTTT' in guide_sequence:
            return "Medium"
        else:
            return "Low"

class InventoryManager:
    """Laboratory inventory management system"""
    
    def __init__(self, db: BioDatabase):
        self.db = db
    
    def check_stock(self, item_name: str) -> Dict:
        """Check inventory stock levels"""
        cursor = self.db.conn.cursor()
        cursor.execute(
            "SELECT quantity, reorder_level FROM inventory WHERE item_name = ?",
            (item_name,)
        )
        result = cursor.fetchone()
        
        if result:
            quantity, reorder_level = result
            status = "Low" if quantity <= reorder_level else "Adequate"
            return {
                'item': item_name,
                'quantity': quantity,
                'reorder_level': reorder_level,
                'status': status,
                'needs_reorder': quantity <= reorder_level
            }
        return None
    
    def add_item(self, item_data: Dict):
        """Add new item to inventory"""
        cursor = self.db.conn.cursor()
        cursor.execute('''
            INSERT INTO inventory (id, item_name, catalog_number, supplier, 
                                  quantity, unit, storage_temp, expiration_date, 
                                  location, reorder_level)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            str(uuid.uuid4()),
            item_data['item_name'],
            item_data.get('catalog_number'),
            item_data.get('supplier'),
            item_data['quantity'],
            item_data['unit'],
            item_data.get('storage_temp'),
            item_data.get('expiration_date'),
            item_data.get('location'),
            item_data.get('reorder_level', 10)
        ))
        self.db.conn.commit()

# ==================== DASHBOARD APPLICATION ====================

class BioLabDashboard:
    """Web-based dashboard for Biolab-OS"""
    
    def __init__(self):
        self.app = dash.Dash(
            __name__,
            external_stylesheets=[dbc.themes.DARKLY],
            suppress_callback_exceptions=True
        )
        self.db = BioDatabase()
        self.setup_layout()
        self.setup_callbacks()
    
    def setup_layout(self):
        """Setup the dashboard layout"""
        self.app.layout = dbc.Container([
            # Navigation Bar
            dbc.NavbarSimple(
                brand="BIOLAB-OS",
                color="primary",
                dark=True,
                children=[
                    dbc.NavItem(dbc.NavLink("Dashboard", href="/")),
                    dbc.NavItem(dbc.NavLink("Virology", href="/virology")),
                    dbc.NavItem(dbc.NavLink("Genetics", href="/genetics")),
                    dbc.NavItem(dbc.NavLink("Inventory", href="/inventory")),
                    dbc.NavItem(dbc.NavLink("Experiments", href="/experiments")),
                ],
                fluid=True,
            ),
            
            # Content
            dcc.Location(id='url', refresh=False),
            html.Div(id='page-content', style={'padding': '20px'}),
            
            # Hidden divs for data storage
            dcc.Store(id='session-data'),
            dcc.Store(id='analysis-results'),
        ], fluid=True)
    
    def create_home_page(self):
        """Create home page with overview"""
        return html.Div([
            dbc.Row([
                dbc.Col(html.H1("Biolab Operating System"), width=12),
                dbc.Col(html.P("Integrated platform for virology, genetics, and molecular biology laboratories"), width=12),
            ]),
            
            dbc.Row([
                # System Status
                dbc.Col(dbc.Card([
                    dbc.CardHeader("System Status"),
                    dbc.CardBody([
                        html.H4("🟢 Operational", className="card-title"),
                        html.P("All systems normal"),
                        html.P(f"Database: Connected"),
                        html.P(f"Users: 5 active"),
                    ])
                ]), width=4),
                
                # Recent Experiments
                dbc.Col(dbc.Card([
                    dbc.CardHeader("Recent Experiments"),
                    dbc.CardBody([
                        html.Ul([
                            html.Li("PCR - SARS-CoV-2 detection"),
                            html.Li("CRISPR - Gene knockout"),
                            html.Li("Sequencing - Viral genome"),
                            html.Li("ELISA - Antibody detection"),
                        ])
                    ])
                ]), width=4),
                
                # Inventory Status
                dbc.Col(dbc.Card([
                    dbc.CardHeader("Inventory Alert"),
                    dbc.CardBody([
                        html.H4("⚠️ 3 items low", className="card-title text-warning"),
                        html.P("Taq Polymerase: 2 vials"),
                        html.P("Primers: 1 set"),
                        html.P("RNase-free tips: 5 boxes"),
                    ])
                ]), width=4),
            ]),
            
            dbc.Row([
                dbc.Col(html.H3("Quick Actions"), width=12),
                dbc.Col(dbc.Button("New Experiment", color="primary", size="lg", className="m-2"), width="auto"),
                dbc.Col(dbc.Button("Add Sample", color="success", size="lg", className="m-2"), width="auto"),
                dbc.Col(dbc.Button("Order Supplies", color="warning", size="lg", className="m-2"), width="auto"),
                dbc.Col(dbc.Button("Generate Report", color="info", size="lg", className="m-2"), width="auto"),
            ]),
            
            # Charts
            dbc.Row([
                dbc.Col(dcc.Graph(id='experiment-stats'), width=6),
                dbc.Col(dcc.Graph(id='inventory-chart'), width=6),
            ]),
        ])
    
    def create_virology_page(self):
        """Create virology analysis page"""
        return html.Div([
            html.H1("Virology Module"),
            
            dbc.Tabs([
                dbc.Tab(label="Virus Database", children=[
                    dbc.Card([
                        dbc.CardBody([
                            html.H4("Virus Catalog"),
                            dbc.Table([
                                html.Thead(html.Tr([
                                    html.Th("Name"),
                                    html.Th("Family"),
                                    html.Th("Genome Type"),
                                    html.Th("Risk Level"),
                                    html.Th("Actions"),
                                ])),
                                html.Tbody([
                                    html.Tr([
                                        html.Td("SARS-CoV-2"),
                                        html.Td("Coronaviridae"),
                                        html.Td("ssRNA(+)"),
                                        html.Td("BSL-3"),
                                        html.Td(dbc.Button("Analyze", size="sm")),
                                    ]),
                                    html.Tr([
                                        html.Td("Influenza A"),
                                        html.Td("Orthomyxoviridae"),
                                        html.Td("ssRNA(-)"),
                                        html.Td("BSL-2"),
                                        html.Td(dbc.Button("Analyze", size="sm")),
                                    ]),
                                ])
                            ])
                        ])
                    ])
                ]),
                
                dbc.Tab(label="Genome Analysis", children=[
                    dbc.Row([
                        dbc.Col([
                            dbc.Textarea(
                                id="virus-sequence-input",
                                placeholder="Paste viral genome sequence here...",
                                style={'width': '100%', 'height': 200},
                            ),
                            dbc.Button("Analyze Genome", id="analyze-virus-btn", className="mt-2"),
                        ], width=6),
                        
                        dbc.Col([
                            html.Div(id="virus-analysis-results"),
                        ], width=6),
                    ]),
                    
                    html.Div(id="virus-charts"),
                ]),
                
                dbbc.Tab(label="PCR Design", children=[
                    dbc.Row([
                        dbc.Col([
                            html.H5("Target Sequence"),
                            dbc.Textarea(id="pcr-target-sequence", style={'width': '100%', 'height': 100}),
                            html.H5("Forward Primer"),
                            dbc.Input(id="forward-primer", type="text"),
                            html.H5("Reverse Primer"),
                            dbc.Input(id="reverse-primer", type="text"),
                            dbc.Button("Simulate PCR", id="simulate-pcr-btn", className="mt-2"),
                        ], width=6),
                        
                        dbc.Col([
                            html.Div(id="pcr-results"),
                        ], width=6),
                    ]),
                ]),
            ]),
        ])
    
    def create_genetics_page(self):
        """Create genetics analysis page"""
        return html.Div([
            html.H1("Genetics Module"),
            
            dbc.Tabs([
                dbc.Tab(label="DNA Analysis", children=[
                    dbc.Row([
                        dbc.Col([
                            dbc.Textarea(
                                id="dna-sequence-input",
                                placeholder="Paste DNA sequence here...",
                                style={'width': '100%', 'height': 200},
                            ),
                            dbc.Button("Analyze DNA", id="analyze-dna-btn", className="mt-2"),
                        ], width=6),
                        
                        dbc.Col([
                            html.Div(id="dna-analysis-results"),
                            html.Div(id="restriction-sites"),
                        ], width=6),
                    ]),
                ]),
                
                dbc.Tab(label="CRISPR Design", children=[
                    dbc.Row([
                        dbc.Col([
                            dbc.Textarea(
                                id="crispr-target-sequence",
                                placeholder="Paste target DNA sequence for CRISPR...",
                                style={'width': '100%', 'height': 150},
                            ),
                            dbc.Select(
                                id="cas-type",
                                options=[
                                    {'label': 'Cas9', 'value': 'cas9'},
                                    {'label': 'Cas12a', 'value': 'cas12a'},
                                    {'label': 'Cas13', 'value': 'cas13'},
                                ],
                                value='cas9'
                            ),
                            dbc.Button("Design gRNA", id="design-crispr-btn", className="mt-2"),
                        ], width=6),
                        
                        dbc.Col([
                            html.Div(id="crispr-design-results"),
                        ], width=6),
                    ]),
                ]),
                
                dbc.Tab(label="Sequence Alignment", children=[
                    dbc.Row([
                        dbc.Col([
                            html.H5("Sequence 1"),
                            dbc.Textarea(id="align-seq1", style={'width': '100%', 'height': 100}),
                            html.H5("Sequence 2"),
                            dbc.Textarea(id="align-seq2", style={'width': '100%', 'height': 100}),
                            dbc.Button("Align Sequences", id="align-btn", className="mt-2"),
                        ], width=6),
                        
                        dbc.Col([
                            html.Div(id="alignment-results"),
                            html.Div(id="alignment-visualization"),
                        ], width=6),
                    ]),
                ]),
            ]),
        ])
    
    def setup_callbacks(self):
        """Setup Dash callbacks"""
        
        @self.app.callback(
            Output('page-content', 'children'),
            Input('url', 'pathname')
        )
        def display_page(pathname):
            if pathname == '/virology':
                return self.create_virology_page()
            elif pathname == '/genetics':
                return self.create_genetics_page()
            elif pathname == '/inventory':
                return self.create_inventory_page()
            elif pathname == '/experiments':
                return self.create_experiments_page()
            else:
                return self.create_home_page()
        
        @self.app.callback(
            Output('virus-analysis-results', 'children'),
            Input('analyze-virus-btn', 'n_clicks'),
            State('virus-sequence-input', 'value')
        )
        def analyze_virus(n_clicks, sequence):
            if n_clicks and sequence:
                analyzer = VirusAnalyzer()
                chars = analyzer.calculate_genome_characteristics(sequence)
                orfs = analyzer.predict_open_reading_frames(sequence)
                
                return html.Div([
                    html.H5("Genome Characteristics"),
                    dbc.Table([
                        html.Tr([html.Th("Length"), html.Td(f"{chars['length']} bp")]),
                        html.Tr([html.Th("GC Content"), html.Td(f"{chars['gc_content']:.2f}%")]),
                        html.Tr([html.Th("AT Content"), html.Td(f"{chars['at_content']:.2f}%")]),
                        html.Tr([html.Th("GC Skew"), html.Td(f"{chars['gc_skew']:.3f}")]),
                    ]),
                    html.H5(f"ORFs Found: {len(orfs)}"),
                    html.Div([
                        html.P(f"ORF {i+1}: Position {orf['start']}-{orf['stop']}, Length: {orf['length']} bp")
                        for i, orf in enumerate(orfs[:5])
                    ])
                ])
            return "Enter a sequence and click analyze"
        
        @self.app.callback(
            Output('pcr-results', 'children'),
            Input('simulate-pcr-btn', 'n_clicks'),
            State('pcr-target-sequence', 'value'),
            State('forward-primer', 'value'),
            State('reverse-primer', 'value')
        )
        def simulate_pcr(n_clicks, target, forward, reverse):
            if n_clicks and target and forward and reverse:
                simulator = PCRSimulator()
                result = simulator.simulate_pcr(target, forward, reverse)
                
                if result['success']:
                    return html.Div([
                        html.H5("PCR Simulation Results"),
                        dbc.Table([
                            html.Tr([html.Th("Product Size"), html.Td(f"{result['product_size']} bp")]),
                            html.Tr([html.Th("Forward Binding"), html.Td(f"Position {result['forward_binding_position']}")]),
                            html.Tr([html.Th("Reverse Binding"), html.Td(f"Position {result['reverse_binding_position']}")]),
                            html.Tr([html.Th("Estimated Yield"), html.Td(f"{result['estimated_product_yield']:.2e} copies")]),
                            html.Tr([html.Th("Efficiency"), html.Td(f"{result['efficiency']*100:.1f}%")]),
                        ])
                    ])
                else:
                    return html.Div([
                        html.H5("PCR Simulation Failed"),
                        html.P(result['error'])
                    ])
            return "Enter all parameters and click simulate"
        
        @self.app.callback(
            Output('crispr-design-results', 'children'),
            Input('design-crispr-btn', 'n_clicks'),
            State('crispr-target-sequence', 'value'),
            State('cas-type', 'value')
        )
        def design_crispr(n_clicks, target_sequence, cas_type):
            if n_clicks and target_sequence:
                designer = CRISPRDesigner()
                guides = designer.design_guide_rna(target_sequence)
                
                table_rows = []
                for i, guide in enumerate(guides[:5]):
                    table_rows.append(html.Tr([
                        html.Td(i+1),
                        html.Td(guide['guide_sequence']),
                        html.Td(guide['pam_sequence']),
                        html.Td(f"{guide['gc_content']:.1f}%"),
                        html.Td(f"{guide['specificity_score']:.1f}"),
                        html.Td(guide['off_target_risk']),
                    ]))
                
                return html.Div([
                    html.H5(f"Top 5 gRNA Designs for {cas_type.upper()}"),
                    dbc.Table([
                        html.Thead(html.Tr([
                            html.Th("#"),
                            html.Th("Guide Sequence (20nt)"),
                            html.Th("PAM"),
                            html.Th("GC%"),
                            html.Th("Specificity"),
                            html.Th("Off-target Risk"),
                        ])),
                        html.Tbody(table_rows)
                    ])
                ])
            return "Enter target sequence and click design"

# ==================== API MODULE ====================

from flask import Flask, jsonify, request
from flask_restful import Api, Resource
from flask_cors import CORS

class BioLabAPI:
    """REST API for BioLab-OS"""
    
    def __init__(self):
        self.flask_app = Flask(__name__)
        self.api = Api(self.flask_app)
        CORS(self.flask_app)
        self.db = BioDatabase()
        
        # Register API endpoints
        self.api.add_resource(SamplesAPI, '/api/samples')
        self.api.add_resource(VirusesAPI, '/api/viruses')
        self.api.add_resource(SequencesAPI, '/api/sequences')
        self.api.add_resource(PCRAPI, '/api/pcr')
        self.api.add_resource(InventoryAPI, '/api/inventory')
        self.api.add_resource(AnalysisAPI, '/api/analyze')
        
        @self.flask_app.route('/api/health')
        def health_check():
            return jsonify({'status': 'healthy', 'timestamp': datetime.datetime.now().isoformat()})

class SamplesAPI(Resource):
    """API for sample management"""
    
    def get(self):
        cursor = BioDatabase().conn.cursor()
        cursor.execute("SELECT * FROM samples ORDER BY created_at DESC LIMIT 100")
        samples = cursor.fetchall()
        return jsonify([dict(zip([column[0] for column in cursor.description], sample)) for sample in samples])
    
    def post(self):
        data = request.get_json()
        sample_id = str(uuid.uuid4())
        
        cursor = BioDatabase().conn.cursor()
        cursor.execute('''
            INSERT INTO samples (id, sample_type, source, collection_date, 
                               location, storage_conditions, volume_ul, 
                               concentration_ng_ul, purity_260_280, created_by)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            sample_id,
            data['sample_type'],
            data['source'],
            data['collection_date'],
            data.get('location'),
            data.get('storage_conditions'),
            data.get('volume_ul'),
            data.get('concentration_ng_ul'),
            data.get('purity_260_280'),
            data.get('created_by', 1)
        ))
        BioDatabase().conn.commit()
        
        return jsonify({'id': sample_id, 'message': 'Sample created'})

# ==================== SECURITY MODULE ====================

class BioSecurity:
    """Biosafety and cybersecurity module"""
    
    @staticmethod
    def validate_sequence(sequence: str) -> bool:
        """Validate DNA/RNA sequence"""
        valid_nucleotides = set('ACGTUacgtuNn')
        return all(nuc in valid_nucleotides for nuc in sequence)
    
    @staticmethod
    def check_dangerous_sequences(sequence: str) -> List[str]:
        """Check for dangerous sequences (toxins, pathogens)"""
        dangerous_motifs = [
            'ATG',  # Start codon (basic check)
            # Add actual dangerous sequence motifs here
            # This would connect to databases like GenBank
        ]
        
        warnings = []
        for motif in dangerous_motifs:
            if motif in sequence:
                warnings.append(f"Found motif: {motif}")
        
        return warnings
    
    @staticmethod
    def assess_risk_level(organism: str, sequence: str) -> str:
        """Assess biosafety risk level"""
        high_risk_organisms = ['SARS-CoV-2', 'Ebola', 'Smallpox', 'Anthrax']
        
        if organism in high_risk_organisms:
            return "BSL-3/4 Required"
        
        # Additional checks
        if len(sequence) > 10000:  # Large genomes might be complex
            return "BSL-2 Recommended"
        
        return "BSL-1/2"

# ==================== MAIN APPLICATION ====================

class BioLabOS:
    """Main Biolab Operating System"""
    
    def __init__(self):
        self.db = BioDatabase()
        self.inventory = InventoryManager(self.db)
        self.security = BioSecurity()
        
        # Initialize modules
        self.virus_analyzer = VirusAnalyzer()
        self.dna_analyzer = DNASequenceAnalyzer()
        self.pcr_simulator = PCRSimulator()
        self.crispr_designer = CRISPRDesigner()
        
        print("""
        ╔══════════════════════════════════════════════════════╗
        ║           BIOLAB-OS v1.0 - Initializing             ║
        ║      Biological Laboratory Operating System         ║
        ╚══════════════════════════════════════════════════════╝
        """)
        print("Modules loaded:")
        print("  ✓ Database System")
        print("  ✓ Virology Analyzer")
        print("  ✓ Genetics Toolkit")
        print("  ✓ PCR Simulator")
        print("  ✓ CRISPR Designer")
        print("  ✓ Inventory Manager")
        print("  ✓ Security Module")
    
    def run_dashboard(self):
        """Start the web dashboard"""
        print("\nStarting BioLab Dashboard on http://localhost:8050")
        dashboard = BioLabDashboard()
        dashboard.app.run_server(debug=True, port=8050)
    
    def run_api(self):
        """Start the REST API"""
        print("\nStarting BioLab API on http://localhost:5000")
        api = BioLabAPI()
        api.flask_app.run(debug=True, port=5000)
    
    def command_line_interface(self):
        """Command line interface for BioLab-OS"""
        print("\nBioLab-OS Command Line Interface")
        print("Type 'help' for commands or 'exit' to quit")
        
        while True:
            command = input("\nBIOLAB> ").strip().lower()
            
            if command == 'exit':
                break
            elif command == 'help':
                self.show_help()
            elif command == 'analyze':
                self.cli_analyze()
            elif command == 'inventory':
                self.cli_inventory()
            elif command == 'sample':
                self.cli_add_sample()
            else:
                print("Unknown command. Type 'help' for available commands.")
    
    def show_help(self):
        """Show available CLI commands"""
        print("\nAvailable commands:")
        print("  analyze    - Analyze DNA/RNA sequences")
        print("  inventory  - Check/manage inventory")
        print("  sample     - Add new biological sample")
        print("  pcr        - Design PCR experiments")
        print("  crispr     - Design CRISPR experiments")
        print("  virus      - Analyze viral genomes")
        print("  dashboard  - Start web dashboard")
        print("  api        - Start REST API")
        print("  help       - Show this help")
        print("  exit       - Exit BioLab-OS")

# ==================== EXAMPLE USAGE ====================

if __name__ == "__main__":
    # Initialize the system
    biolab = BioLabOS()
    
    # Example 1: Analyze a viral genome
    print("\n=== Example: Viral Genome Analysis ===")
    test_sequence = "ATGCGATACGTTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    analysis = biolab.virus_analyzer.calculate_genome_characteristics(test_sequence)
    print(f"Genome Length: {analysis['length']} bp")
    print(f"GC Content: {analysis['gc_content']:.2f}%")
    
    # Example 2: Design CRISPR guides
    print("\n=== Example: CRISPR Guide Design ===")
    target_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    guides = biolab.crispr_designer.design_guide_rna(target_seq)
    print(f"Designed {len(guides)} guide RNAs")
    for i, guide in enumerate(guides[:3]):
        print(f"Guide {i+1}: {guide['guide_sequence']} (Score: {guide['specificity_score']:.1f})")
    
    # Example 3: Check inventory
    print("\n=== Example: Inventory Check ===")
    # This would check actual database
    
    # Choose interface
    print("\nSelect interface:")
    print("1. Web Dashboard")
    print("2. REST API")
    print("3. Command Line")
    
    choice = input("Enter choice (1-3): ").strip()
    
    if choice == "1":
        biolab.run_dashboard()
    elif choice == "2":
        biolab.run_api()
    elif choice == "3":
        biolab.command_line_interface()
