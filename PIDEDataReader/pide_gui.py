#!/usr/bin/env python3
"""
PIDE Database Explorer GUI (Tkinter version)

A graphical interface for exploring cell survival data from the 
GSI Particle Irradiation Data Ensemble (PIDE) database.

Features:
- Select cell type and radiation type
- Display raw experimental data as scatter points
- Display LQ model curves
- Show experiment details and LQ parameters

Requirements:
- tkinter (usually comes with Python)
- matplotlib
- pandas
- numpy

Usage:
    python3 pide_gui.py [path_to_PIDE3.4_folder]
"""

import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path

import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure


class PIDEDataManager:
    """Manages loading and querying PIDE data"""
    
    def __init__(self, pide_dir):
        self.pide_dir = Path(pide_dir)
        self.main_df = None
        self.ion_raw_data = {}
        self.photon_raw_data = {}
        self.loaded = False
        
    def load_data(self):
        """Load all PIDE data files"""
        try:
            # Load main CSV database
            csv_file = self.pide_dir / "PIDE3.4.csv"
            if not csv_file.exists():
                raise FileNotFoundError(f"CSV file not found: {csv_file}\nPlease run convert_pide_xlsx.py first.")
            
            self.main_df = pd.read_csv(csv_file)
            
            # Load ion raw data
            ion_raw_file = self.pide_dir / "PIDE3.4_IonRawData.dat"
            if ion_raw_file.exists():
                self._load_ion_raw_data(ion_raw_file)
            
            # Load photon raw data
            photon_raw_file = self.pide_dir / "PIDE3.4_PhotonRawData.dat"
            if photon_raw_file.exists():
                self._load_photon_raw_data(photon_raw_file)
            
            self.loaded = True
            return True
            
        except Exception as e:
            print(f"Error loading data: {e}")
            return False
    
    def _load_ion_raw_data(self, filename):
        """Parse ion raw data file"""
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip().replace('\r', '')
                parts = line.split('\t')
                if len(parts) < 6:
                    continue
                
                exp_id = int(parts[0])
                
                # Check for N/A
                if parts[5] == 'N/A':
                    continue
                
                # Parse dose-survival pairs (starting from column 5)
                doses = []
                sfs = []
                values = parts[5:]
                for i in range(0, len(values)-1, 2):
                    try:
                        doses.append(float(values[i]))
                        sfs.append(float(values[i+1]))
                    except (ValueError, IndexError):
                        break
                
                if len(doses) > 1:
                    self.ion_raw_data[exp_id] = {
                        'doses': np.array(doses),
                        'sf': np.array(sfs)
                    }
    
    def _load_photon_raw_data(self, filename):
        """Parse photon raw data file"""
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip().replace('\r', '')
                parts = line.split('\t')
                if len(parts) < 5:
                    continue
                
                exp_id = int(parts[0])
                pub_id = int(parts[1])
                photon_exp = int(parts[3])
                
                # Check for N/A
                if parts[4] == 'N/A':
                    continue
                
                # Parse dose-survival pairs
                doses = []
                sfs = []
                values = parts[4:]
                for i in range(0, len(values)-1, 2):
                    try:
                        doses.append(float(values[i]))
                        sfs.append(float(values[i+1]))
                    except (ValueError, IndexError):
                        break
                
                if len(doses) > 1:
                    key = (pub_id, photon_exp)
                    self.photon_raw_data[key] = {
                        'doses': np.array(doses),
                        'sf': np.array(sfs)
                    }
    
    def get_unique_cell_lines(self):
        """Get list of unique cell lines"""
        if self.main_df is None:
            return []
        return sorted(self.main_df['Cells'].unique().tolist())
    
    def get_ions_for_cell(self, cell_line):
        """Get list of ions available for a specific cell line"""
        if self.main_df is None:
            return []
        cell_data = self.main_df[self.main_df['Cells'] == cell_line]
        return sorted(cell_data['Ion'].unique().tolist())
    
    def get_experiments(self, cell_line, ion):
        """Get all experiments for a cell line and ion combination"""
        if self.main_df is None:
            return pd.DataFrame()
        return self.main_df[(self.main_df['Cells'] == cell_line) & 
                           (self.main_df['Ion'] == ion)].copy()
    
    def get_unique_ions(self):
        """Get list of unique ion/radiation types"""
        if self.main_df is None:
            return []
        return sorted(self.main_df['Ion'].unique().tolist())
    
    def get_unique_photon_types(self):
        """Get list of unique photon radiation types"""
        if self.main_df is None:
            return []
        photon_types = self.main_df['PhotonRadiation'].dropna().unique()
        return sorted([p for p in photon_types if p and str(p).strip()])
    
    def get_lets_for_ion(self, ion):
        """Get list of unique LET values for a specific ion"""
        if self.main_df is None:
            return []
        ion_data = self.main_df[self.main_df['Ion'] == ion]
        let_values = ion_data['LET'].dropna().unique()
        return sorted(let_values.tolist())
    
    def get_experiments_by_ion_let(self, ion, let=None, cell_filter=None):
        """Get experiments for a specific ion and optionally LET value"""
        if self.main_df is None:
            return pd.DataFrame()
        
        mask = self.main_df['Ion'] == ion
        
        if let is not None:
            # Allow some tolerance for LET matching (±0.5 keV/μm)
            mask = mask & (abs(self.main_df['LET'] - let) < 0.5)
        
        if cell_filter and cell_filter != 'All':
            mask = mask & (self.main_df['Cells'] == cell_filter)
        
        return self.main_df[mask].copy()
    
    def get_cells_for_ion(self, ion):
        """Get list of cell lines that have data for a specific ion"""
        if self.main_df is None:
            return []
        ion_data = self.main_df[self.main_df['Ion'] == ion]
        return sorted(ion_data['Cells'].unique().tolist())
    
    def get_cells_for_ion_let(self, ion, let):
        """Get list of cell lines that have data for a specific ion and LET"""
        if self.main_df is None:
            return []
        mask = (self.main_df['Ion'] == ion) & (abs(self.main_df['LET'] - let) < 0.5)
        return sorted(self.main_df[mask]['Cells'].unique().tolist())
    
    def get_experiments_by_photon(self, photon_type, cell_filter=None):
        """Get experiments for a specific photon radiation type"""
        if self.main_df is None:
            return pd.DataFrame()
        
        mask = self.main_df['PhotonRadiation'] == photon_type
        
        if cell_filter and cell_filter != 'All':
            mask = mask & (self.main_df['Cells'] == cell_filter)
        
        return self.main_df[mask].copy()
    
    def get_cells_for_photon(self, photon_type):
        """Get list of cell lines that have data for a specific photon type"""
        if self.main_df is None:
            return []
        photon_data = self.main_df[self.main_df['PhotonRadiation'] == photon_type]
        return sorted(photon_data['Cells'].unique().tolist())
    
    def get_ion_raw_data(self, exp_id):
        """Get raw survival data for an experiment"""
        return self.ion_raw_data.get(exp_id, None)
    
    def get_photon_raw_data(self, pub_id, photon_exp):
        """Get photon reference raw data"""
        return self.photon_raw_data.get((pub_id, photon_exp), None)


class PIDEExplorerGUI:
    """Main GUI class for PIDE database exploration"""
    
    def __init__(self, root, pide_dir):
        self.root = root
        self.root.title('PIDE Database Explorer - Cell Survival Curves')
        self.root.geometry('1500x950')
        
        self.pide_dir = pide_dir
        self.data_manager = PIDEDataManager(pide_dir)
        self.current_experiments = pd.DataFrame()
        # Track plotted elements for highlighting
        self.plot_elements = {}
        self.selected_exp_id = None
        self.highlight_annotation = None
        self.temp_highlight_elements = []  # Temporary elements created during highlighting
        # Track fitted parameters for annotation
        self.params_fitted = False
        self.fitted_alpha = None
        self.fitted_beta = None
        
        # Variables for checkboxes
        self.show_raw_var = tk.BooleanVar(value=True)
        self.show_lq_var = tk.BooleanVar(value=True)
        self.show_photon_var = tk.BooleanVar(value=False)
        
        # Selection mode variable
        self.selection_mode = tk.StringVar(value='cell_ion')  # 'cell_ion' or 'radiation_let'
        
        self.create_widgets()
        self.load_data()
        
    def create_widgets(self):
        """Create all GUI widgets"""
        # Main container
        main_frame = ttk.Frame(self.root, padding="5")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Left panel - Controls
        left_frame = ttk.Frame(main_frame, width=300)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_frame.pack_propagate(False)
        
        # Right panel - Plot and table
        right_frame = ttk.Frame(main_frame)
        right_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # === Left Panel Contents ===
        
        # Title
        title_label = ttk.Label(left_frame, text='PIDE Database Explorer', 
                               font=('Arial', 14, 'bold'))
        title_label.pack(pady=10)
        
        # Selection Mode
        mode_frame = ttk.LabelFrame(left_frame, text='Selection Mode', padding=10)
        mode_frame.pack(fill=tk.X, pady=5)
        
        ttk.Radiobutton(mode_frame, text='By Cell Line → Ion', 
                       variable=self.selection_mode, value='cell_ion',
                       command=self.on_mode_changed).pack(anchor=tk.W)
        ttk.Radiobutton(mode_frame, text='By Ion Type → LET', 
                       variable=self.selection_mode, value='radiation_let',
                       command=self.on_mode_changed).pack(anchor=tk.W)
        ttk.Radiobutton(mode_frame, text='By Photon Type', 
                       variable=self.selection_mode, value='photon',
                       command=self.on_mode_changed).pack(anchor=tk.W)
        
        # Cell type selection (Mode 1)
        self.cell_frame = ttk.LabelFrame(left_frame, text='Cell Type Selection', padding=10)
        self.cell_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(self.cell_frame, text='Cell Line:').pack(anchor=tk.W)
        self.cell_combo = ttk.Combobox(self.cell_frame, state='readonly', width=25)
        self.cell_combo.pack(fill=tk.X, pady=2)
        self.cell_combo.bind('<<ComboboxSelected>>', self.on_cell_changed)
        
        # Ion selection (Mode 1)
        self.ion_frame = ttk.LabelFrame(left_frame, text='Ion Selection', padding=10)
        self.ion_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(self.ion_frame, text='Ion Type:').pack(anchor=tk.W)
        self.ion_combo = ttk.Combobox(self.ion_frame, state='readonly', width=25)
        self.ion_combo.pack(fill=tk.X, pady=2)
        self.ion_combo.bind('<<ComboboxSelected>>', self.on_ion_changed)
        
        # Radiation Type selection (Mode 2)
        self.rad_type_frame = ttk.LabelFrame(left_frame, text='Radiation Type', padding=10)
        
        ttk.Label(self.rad_type_frame, text='Ion/Particle:').pack(anchor=tk.W)
        self.rad_type_combo = ttk.Combobox(self.rad_type_frame, state='readonly', width=25)
        self.rad_type_combo.pack(fill=tk.X, pady=2)
        self.rad_type_combo.bind('<<ComboboxSelected>>', self.on_rad_type_changed)
        
        # LET selection (Mode 2)
        self.let_frame = ttk.LabelFrame(left_frame, text='LET Selection', padding=10)
        
        ttk.Label(self.let_frame, text='LET (keV/μm):').pack(anchor=tk.W)
        self.let_combo = ttk.Combobox(self.let_frame, state='readonly', width=25)
        self.let_combo.pack(fill=tk.X, pady=2)
        self.let_combo.bind('<<ComboboxSelected>>', self.on_let_changed)
        
        # LET range info
        self.let_info_label = ttk.Label(self.let_frame, text='', foreground='blue', wraplength=250)
        self.let_info_label.pack(anchor=tk.W, pady=5)
        
        # Cell filter for Mode 2 (optional)
        self.cell_filter_frame = ttk.LabelFrame(left_frame, text='Cell Filter (Optional)', padding=10)
        
        ttk.Label(self.cell_filter_frame, text='Filter by Cell:').pack(anchor=tk.W)
        self.cell_filter_combo = ttk.Combobox(self.cell_filter_frame, state='readonly', width=25)
        self.cell_filter_combo.pack(fill=tk.X, pady=2)
        self.cell_filter_combo.bind('<<ComboboxSelected>>', self.on_cell_filter_changed)
        
        # Photon Type selection (Mode 3)
        self.photon_type_frame = ttk.LabelFrame(left_frame, text='Photon Radiation Type', padding=10)
        
        ttk.Label(self.photon_type_frame, text='Photon Source:').pack(anchor=tk.W)
        self.photon_type_combo = ttk.Combobox(self.photon_type_frame, state='readonly', width=25)
        self.photon_type_combo.pack(fill=tk.X, pady=2)
        self.photon_type_combo.bind('<<ComboboxSelected>>', self.on_photon_type_changed)
        
        # Photon info label
        self.photon_info_label = ttk.Label(self.photon_type_frame, text='', foreground='blue', wraplength=250)
        self.photon_info_label.pack(anchor=tk.W, pady=5)
        
        # Cell filter for Mode 3 (optional)
        self.photon_cell_filter_frame = ttk.LabelFrame(left_frame, text='Cell Filter (Optional)', padding=10)
        
        ttk.Label(self.photon_cell_filter_frame, text='Filter by Cell:').pack(anchor=tk.W)
        self.photon_cell_filter_combo = ttk.Combobox(self.photon_cell_filter_frame, state='readonly', width=25)
        self.photon_cell_filter_combo.pack(fill=tk.X, pady=2)
        self.photon_cell_filter_combo.bind('<<ComboboxSelected>>', self.on_photon_cell_filter_changed)
        
        # Display options
        display_frame = ttk.LabelFrame(left_frame, text='Display Options', padding=10)
        display_frame.pack(fill=tk.X, pady=5)
        
        ttk.Checkbutton(display_frame, text='Show Raw Data (dots)', 
                       variable=self.show_raw_var, 
                       command=self.update_plot).pack(anchor=tk.W)
        ttk.Checkbutton(display_frame, text='Show LQ Curves (lines)', 
                       variable=self.show_lq_var,
                       command=self.update_plot).pack(anchor=tk.W)
        ttk.Checkbutton(display_frame, text='Show Photon Reference', 
                       variable=self.show_photon_var,
                       command=self.update_plot).pack(anchor=tk.W)
        
        # Statistics
        stats_frame = ttk.LabelFrame(left_frame, text='Statistics', padding=10)
        stats_frame.pack(fill=tk.X, pady=5)
        
        self.stats_label = ttk.Label(stats_frame, text='No data loaded', wraplength=250)
        self.stats_label.pack(anchor=tk.W)
        
        # Buttons
        btn_frame = ttk.Frame(left_frame)
        btn_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(btn_frame, text='Refresh', command=self.update_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text='Clear', command=self.clear_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text='Save Plot', command=self.save_plot).pack(side=tk.LEFT, padx=5)
        
        # Info
        info_label = ttk.Label(left_frame, text='PIDE Database v3.4\nGSI Biophysics',
                              foreground='gray')
        info_label.pack(side=tk.BOTTOM, pady=10)
        
        # === Right Panel Contents ===
        
        # Plot frame
        plot_frame = ttk.Frame(right_frame)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.setup_plot()
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        
        # Toolbar
        toolbar_frame = ttk.Frame(plot_frame)
        toolbar_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()
        
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Table frame
        table_frame = ttk.LabelFrame(right_frame, text='Experiment Details', padding=5)
        table_frame.pack(fill=tk.X, pady=5)
        
        # Create treeview for experiment table
        columns = ('ExpID', 'Cell', 'Publication', 'LET', 'Energy', 'Alpha', 'Beta', 'D10', 'PhotonRef', 'HasRaw')
        self.tree = ttk.Treeview(table_frame, columns=columns, show='headings', height=8)
        
        # Column headings
        self.tree.heading('ExpID', text='ExpID')
        self.tree.heading('Cell', text='Cell Line')
        self.tree.heading('Publication', text='Publication')
        self.tree.heading('LET', text='LET (keV/μm)')
        self.tree.heading('Energy', text='Energy (MeV/u)')
        self.tree.heading('Alpha', text='α (Gy⁻¹)')
        self.tree.heading('Beta', text='β (Gy⁻²)')
        self.tree.heading('D10', text='D10 (Gy)')
        self.tree.heading('PhotonRef', text='Photon Ref')
        self.tree.heading('HasRaw', text='Raw Data')
        
        # Column widths
        self.tree.column('ExpID', width=50)
        self.tree.column('Cell', width=70)
        self.tree.column('Publication', width=80)
        self.tree.column('LET', width=70)
        self.tree.column('Energy', width=70)
        self.tree.column('Alpha', width=60)
        self.tree.column('Beta', width=60)
        self.tree.column('D10', width=60)
        self.tree.column('PhotonRef', width=70)
        self.tree.column('HasRaw', width=55)
        
        # Bind selection event
        self.tree.bind('<<TreeviewSelect>>', self.on_table_select)
        
        # Scrollbar
        scrollbar = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Status bar
        self.status_var = tk.StringVar(value='Ready')
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X, side=tk.BOTTOM)
        
        # Keyboard shortcuts
        self.root.bind('<Escape>', lambda e: self.clear_highlight())
    
    def setup_plot(self):
        """Initialize plot settings"""
        self.ax.set_xlabel('Dose (Gy)', fontsize=12)
        self.ax.set_ylabel('Survival Fraction', fontsize=12)
        self.ax.set_yscale('log')
        self.ax.set_xlim(0, 12)
        self.ax.set_ylim(1e-4, 1.5)
        self.ax.grid(True, alpha=0.3)
        self.ax.set_title('Cell Survival Curve', fontsize=14)
        self.fig.tight_layout()
    
    def load_data(self):
        """Load PIDE data"""
        self.status_var.set('Loading PIDE data...')
        self.root.update()
        
        if self.data_manager.load_data():
            # Populate cell combo (Mode 1)
            cell_lines = self.data_manager.get_unique_cell_lines()
            self.cell_combo['values'] = cell_lines
            
            # Populate radiation type combo (Mode 2)
            rad_types = self.data_manager.get_unique_ions()
            self.rad_type_combo['values'] = rad_types
            
            # Populate photon type combo (Mode 3)
            photon_types = self.data_manager.get_unique_photon_types()
            self.photon_type_combo['values'] = photon_types
            
            # Select V79 by default if available (Mode 1)
            if 'V79' in cell_lines:
                self.cell_combo.set('V79')
                self.on_cell_changed(None)
            elif cell_lines:
                self.cell_combo.set(cell_lines[0])
                self.on_cell_changed(None)
            
            # Initialize mode display
            self.on_mode_changed()
            
            self.status_var.set(f'Loaded {len(self.data_manager.main_df)} experiments')
        else:
            messagebox.showerror('Error', 
                               f'Failed to load PIDE data from:\n{self.pide_dir}\n\n'
                               'Please make sure the PIDE3.4 folder exists and contains:\n'
                               '- PIDE3.4.csv (run convert_pide_xlsx.py first)\n'
                               '- PIDE3.4_IonRawData.dat\n'
                               '- PIDE3.4_PhotonRawData.dat')
            self.status_var.set('Failed to load data')
    
    def on_cell_changed(self, event):
        """Handle cell line selection change"""
        cell_line = self.cell_combo.get()
        if not cell_line:
            return
        
        # Update ion combo
        ions = self.data_manager.get_ions_for_cell(cell_line)
        self.ion_combo['values'] = ions
        
        # Select 12C by default if available
        if '12C' in ions:
            self.ion_combo.set('12C')
        elif ions:
            self.ion_combo.set(ions[0])
        
        self.on_ion_changed(None)
    
    def on_ion_changed(self, event):
        """Handle ion selection change"""
        self.update_plot()
    
    def on_mode_changed(self):
        """Handle selection mode change"""
        mode = self.selection_mode.get()
        
        # Hide all mode-specific widgets first
        self.cell_frame.pack_forget()
        self.ion_frame.pack_forget()
        self.rad_type_frame.pack_forget()
        self.let_frame.pack_forget()
        self.cell_filter_frame.pack_forget()
        self.photon_type_frame.pack_forget()
        self.photon_cell_filter_frame.pack_forget()
        
        # Get the mode frame to pack widgets after
        mode_frame = None
        for child in self.cell_frame.master.winfo_children():
            if isinstance(child, ttk.LabelFrame) and 'Mode' in str(child.cget('text')):
                mode_frame = child
                break
        
        if mode == 'cell_ion':
            # Show Mode 1 widgets
            self.cell_frame.pack(fill=tk.X, pady=5)
            self.ion_frame.pack(fill=tk.X, pady=5)
            
            # Update plot with Mode 1 data
            self.update_plot()
            
        elif mode == 'radiation_let':
            # Show Mode 2 widgets
            self.rad_type_frame.pack(fill=tk.X, pady=5)
            self.let_frame.pack(fill=tk.X, pady=5)
            self.cell_filter_frame.pack(fill=tk.X, pady=5)
            
            # Initialize radiation type if not set
            if not self.rad_type_combo.get() and self.rad_type_combo['values']:
                # Default to 12C if available
                rad_types = list(self.rad_type_combo['values'])
                if '12C' in rad_types:
                    self.rad_type_combo.set('12C')
                else:
                    self.rad_type_combo.set(rad_types[0])
                self.on_rad_type_changed(None)
        
        elif mode == 'photon':
            # Show Mode 3 widgets
            self.photon_type_frame.pack(fill=tk.X, pady=5)
            self.photon_cell_filter_frame.pack(fill=tk.X, pady=5)
            
            # Initialize photon type if not set
            if not self.photon_type_combo.get() and self.photon_type_combo['values']:
                photon_types = list(self.photon_type_combo['values'])
                # Default to 60Co if available
                if '60Co' in photon_types:
                    self.photon_type_combo.set('60Co')
                else:
                    self.photon_type_combo.set(photon_types[0])
                self.on_photon_type_changed(None)
    
    def on_rad_type_changed(self, event):
        """Handle radiation type selection change (Mode 2)"""
        rad_type = self.rad_type_combo.get()
        if not rad_type:
            return
        
        # Update LET combo
        let_values = self.data_manager.get_lets_for_ion(rad_type)
        let_strings = [f'{let:.2f}' for let in let_values]
        self.let_combo['values'] = ['All'] + let_strings
        
        # Update LET range info
        if let_values:
            self.let_info_label.config(
                text=f'LET range: {min(let_values):.1f} - {max(let_values):.1f} keV/μm\n'
                     f'({len(let_values)} unique values)')
        else:
            self.let_info_label.config(text='No LET data available')
        
        # Update cell filter combo
        cells = self.data_manager.get_cells_for_ion(rad_type)
        self.cell_filter_combo['values'] = ['All'] + cells
        self.cell_filter_combo.set('All')
        
        # Select first LET by default
        if let_strings:
            self.let_combo.set('All')
        
        self.on_let_changed(None)
    
    def on_let_changed(self, event):
        """Handle LET selection change (Mode 2)"""
        rad_type = self.rad_type_combo.get()
        let_str = self.let_combo.get()
        
        if not rad_type:
            return
        
        # Update cell filter based on selected LET
        if let_str and let_str != 'All':
            let_val = float(let_str)
            cells = self.data_manager.get_cells_for_ion_let(rad_type, let_val)
            current_filter = self.cell_filter_combo.get()
            self.cell_filter_combo['values'] = ['All'] + cells
            if current_filter not in ['All'] + cells:
                self.cell_filter_combo.set('All')
        
        self.update_plot_mode2()
    
    def on_cell_filter_changed(self, event):
        """Handle cell filter change (Mode 2)"""
        self.update_plot_mode2()
    
    def on_photon_type_changed(self, event):
        """Handle photon type selection change (Mode 3)"""
        photon_type = self.photon_type_combo.get()
        if not photon_type:
            return
        
        # Update cell filter combo
        cells = self.data_manager.get_cells_for_photon(photon_type)
        self.photon_cell_filter_combo['values'] = ['All'] + cells
        self.photon_cell_filter_combo.set('All')
        
        # Update info label
        n_exp = len(self.data_manager.get_experiments_by_photon(photon_type))
        self.photon_info_label.config(
            text=f'{n_exp} experiments with {photon_type}\n'
                 f'{len(cells)} cell lines available')
        
        self.update_plot_mode3()
    
    def on_photon_cell_filter_changed(self, event):
        """Handle cell filter change (Mode 3)"""
        self.update_plot_mode3()

    def on_table_select(self, event):
        """Handle table row selection - highlight corresponding curve"""
        selection = self.tree.selection()
        if not selection:
            self.clear_highlight()
            return

        item = self.tree.item(selection[0])
        if not item.get('values'):
            self.clear_highlight()
            return

        exp_id = int(item['values'][0])
        self.highlight_experiment(exp_id)

    def highlight_experiment(self, exp_id):
        """Highlight the curve for a specific experiment
        
        Always shows both big dots AND thick solid line for the selected experiment,
        creating temporary elements if the original data only had one or the other.
        """
        self.selected_exp_id = exp_id
        
        # Remove any temporary highlight elements from previous selection
        for artist in self.temp_highlight_elements:
            try:
                artist.remove()
            except:
                pass
        self.temp_highlight_elements = []
        
        # Get the selected experiment data for creating temporary elements
        exp_data = self.current_experiments[self.current_experiments['ExpID'] == exp_id]
        exp_row = exp_data.iloc[0] if not exp_data.empty else None

        for eid, elements in self.plot_elements.items():
            is_selected = (eid == exp_id)

            scatter_artist = elements.get('scatter')
            line_artist = elements.get('line')
            
            if scatter_artist is not None:
                scatter_artist.set_alpha(1.0 if is_selected else 0.15)
                scatter_artist.set_sizes([120 if is_selected else 50])
                if is_selected:
                    scatter_artist.set_edgecolors('black')
                    scatter_artist.set_linewidths(2)
                else:
                    scatter_artist.set_edgecolors('none')
                    scatter_artist.set_linewidths(0)

            if line_artist is not None:
                line_artist.set_alpha(1.0 if is_selected else 0.15)
                line_artist.set_linewidth(4 if is_selected else 1.5)
                if is_selected:
                    line_artist.set_linestyle('-')  # Always solid when selected
            
            # For the selected experiment, create missing elements
            if is_selected and exp_row is not None:
                color = scatter_artist.get_facecolors()[0] if scatter_artist is not None else \
                        (line_artist.get_color() if line_artist is not None else 'blue')
                
                # Determine parameter mode based on current selection mode
                current_mode = self.selection_mode.get()
                param_mode = 'photon' if current_mode == 'photon' else 'ion'
                
                # Get LQ parameters and check if complete (both alpha AND beta reported)
                alpha_val, beta_val = self._get_lq_params(exp_row, mode=param_mode)
                # Parameters are complete only if both are reported (beta_val != 0 means it was actually reported)
                params_complete = (alpha_val is not None) and (beta_val is not None and beta_val != 0)
                self.params_fitted = False
                self.fitted_alpha = None
                self.fitted_beta = None
                
                # If parameters incomplete but have raw data, fit the model
                if not params_complete and scatter_artist is not None:
                    offsets = scatter_artist.get_offsets()
                    doses_raw = offsets[:, 0]
                    sf_raw = offsets[:, 1]
                    fitted_alpha, fitted_beta = self._fit_lq_model(doses_raw, sf_raw)
                    if fitted_alpha is not None:
                        alpha_val = fitted_alpha
                        beta_val = fitted_beta
                        self.params_fitted = True
                        self.fitted_alpha = fitted_alpha
                        self.fitted_beta = fitted_beta
                
                # If we fitted parameters, hide existing line and create new fitted line
                if self.params_fitted and alpha_val is not None:
                    # Hide existing line if present (will be restored on clear)
                    if line_artist is not None:
                        line_artist.set_alpha(0)
                    # Create temporary line with fitted parameters
                    doses = np.linspace(0, 12, 100)
                    sf = np.exp(-alpha_val * doses - beta_val * doses**2)
                    temp_line, = self.ax.plot(
                        doses, sf, color=color, linestyle='-',
                        linewidth=4, alpha=1.0, zorder=10
                    )
                    self.temp_highlight_elements.append(temp_line)
                elif line_artist is None and alpha_val is not None:
                    # No existing line, create temporary line from LQ parameters
                    doses = np.linspace(0, 12, 100)
                    sf = np.exp(-alpha_val * doses - beta_val * doses**2)
                    temp_line, = self.ax.plot(
                        doses, sf, color=color, linestyle='-',
                        linewidth=4, alpha=1.0, zorder=10
                    )
                    self.temp_highlight_elements.append(temp_line)
                
                # If no scatter exists but we have line, create temporary dots
                if scatter_artist is None and line_artist is not None:
                    # Get points from the line data
                    line_x = line_artist.get_xdata()
                    line_y = line_artist.get_ydata()
                    # Sample key points (0, 2, 4, 6, 8, 10 Gy)
                    key_doses = [0, 2, 4, 6, 8, 10]
                    key_sf = []
                    for d in key_doses:
                        idx = int(d / 12 * (len(line_x) - 1))
                        key_sf.append(line_y[idx])
                    temp_scatter = self.ax.scatter(
                        key_doses, key_sf,
                        color=color, s=120, alpha=1.0, zorder=10,
                        edgecolors='black', linewidths=2
                    )
                    self.temp_highlight_elements.append(temp_scatter)

        self.add_highlight_annotation(exp_id)
        self.canvas.draw()
    
    def _get_lq_params(self, exp_row, mode='ion'):
        """Get consistent alpha/beta parameters for LQ model.
        
        Args:
            exp_row: DataFrame row with experiment data
            mode: 'ion' or 'photon' - determines which parameter set to use
        
        Returns:
            (alpha, beta) tuple. If beta is N/A, returns 0 for beta.
            Returns (None, None) if no parameters available.
        """
        if mode == 'ion':
            # Use ion parameters (ai_paper, bi_paper)
            if pd.notna(exp_row.get('ai_paper')):
                alpha = exp_row['ai_paper']
                beta = exp_row.get('bi_paper') if pd.notna(exp_row.get('bi_paper')) else 0
                return alpha, beta
            # Fallback to fitted ion parameters
            if pd.notna(exp_row.get('ai_fit')):
                alpha = exp_row['ai_fit']
                beta = exp_row.get('bi_fit') if pd.notna(exp_row.get('bi_fit')) else 0
                return alpha, beta
        else:  # mode == 'photon'
            # Use photon parameters (ax_paper, bx_paper)
            if pd.notna(exp_row.get('ax_paper')):
                alpha = exp_row['ax_paper']
                beta = exp_row.get('bx_paper') if pd.notna(exp_row.get('bx_paper')) else 0
                return alpha, beta
            # Fallback to fitted photon parameters
            if pd.notna(exp_row.get('ax_fit')):
                alpha = exp_row['ax_fit']
                beta = exp_row.get('bx_fit') if pd.notna(exp_row.get('bx_fit')) else 0
                return alpha, beta
        
        return None, None
    
    def _fit_lq_model(self, doses, sf_values):
        """Fit LQ model to raw survival data using numpy linear algebra.
        
        Model: SF = exp(-alpha*D - beta*D^2)
        Taking ln: ln(SF) = -alpha*D - beta*D^2
        This is linear in alpha and beta, solvable via least squares.
        
        Args:
            doses: array of dose values (Gy)
            sf_values: array of survival fraction values
        
        Returns:
            (alpha, beta) tuple, or (None, None) if fitting fails
        """
        try:
            # Filter out any invalid values (D >= 0, SF > 0 for log)
            doses_arr = np.array(doses)
            sf_arr = np.array(sf_values)
            valid_mask = (doses_arr >= 0) & (sf_arr > 0)
            doses_valid = doses_arr[valid_mask]
            sf_valid = sf_arr[valid_mask]
            
            if len(doses_valid) < 2:
                return None, None
            
            # Transform to linear problem: ln(SF) = -alpha*D - beta*D^2
            # Let y = ln(SF), x1 = D, x2 = D^2
            # Then y = -alpha*x1 - beta*x2
            y = np.log(sf_valid)
            X = np.column_stack([doses_valid, doses_valid**2])
            
            # Solve using least squares: X @ params = y
            # params = [−alpha, −beta]
            params, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
            alpha = -params[0]
            beta = -params[1]
            
            # Ensure non-negative values (physically meaningful)
            alpha = max(0, alpha)
            beta = max(0, beta)
            
            return alpha, beta
        except:
            return None, None

    def add_highlight_annotation(self, exp_id):
        """Add text annotation for highlighted experiment"""
        if self.highlight_annotation:
            self.highlight_annotation.remove()
            self.highlight_annotation = None

        exp = self.current_experiments[self.current_experiments['ExpID'] == exp_id]
        if exp.empty:
            return

        exp = exp.iloc[0]
        cell = exp['Cells']
        
        # Check if we used curve-fitted parameters (fitted to raw data during highlight)
        if hasattr(self, 'params_fitted') and self.params_fitted:
            # Use fitted parameters
            alpha = self.fitted_alpha
            beta = self.fitted_beta
            text = f"ExpID: {exp_id}\nCell: {cell}"
            text += f"\nα={alpha:.4f}, β={beta:.5f}"
            text += "\n(curve fitted)"
        else:
            # Track source of parameters
            source = None
            
            # Use ion parameters if available (ai_paper first, then ai_fit)
            if pd.notna(exp.get('ai_paper')):
                alpha = exp['ai_paper']
                beta = exp.get('bi_paper')
                source = 'paper'
            elif pd.notna(exp.get('ai_fit')):
                alpha = exp['ai_fit']
                beta = exp.get('bi_fit')
                source = 'fit'
            # Otherwise try photon parameters (ax_paper first, then ax_fit)
            elif pd.notna(exp.get('ax_paper')):
                alpha = exp['ax_paper']
                beta = exp.get('bx_paper')
                source = 'paper'
            elif pd.notna(exp.get('ax_fit')):
                alpha = exp['ax_fit']
                beta = exp.get('bx_fit')
                source = 'fit'
            else:
                alpha = None
                beta = None
                source = None

            text = f"ExpID: {exp_id}\nCell: {cell}"
            if pd.notna(alpha):
                text += f"\nα={alpha:.4f}"
                if pd.notna(beta):
                    text += f", β={beta:.5f}"
                if source == 'fit':
                    text += "\n(from fit)"

        self.highlight_annotation = self.ax.annotate(
            text,
            xy=(0.98, 0.98),
            xycoords='axes fraction',
            ha='right',
            va='top',
            fontsize=9,
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.9, edgecolor='black')
        )

    def clear_highlight(self):
        """Clear any existing highlight and remove temporary elements"""
        # Remove temporary highlight elements
        for artist in self.temp_highlight_elements:
            try:
                artist.remove()
            except:
                pass
        self.temp_highlight_elements = []
        
        # Reset fitted parameters flag
        self.params_fitted = False
        self.fitted_alpha = None
        self.fitted_beta = None
        
        # Restore original styling for all plot elements
        for elements in self.plot_elements.values():
            scatter_artist = elements.get('scatter')
            if scatter_artist is not None:
                scatter_artist.set_alpha(0.8)
                scatter_artist.set_sizes([50])
                scatter_artist.set_edgecolors('none')
                scatter_artist.set_linewidths(0)

            line_artist = elements.get('line')
            if line_artist is not None:
                line_artist.set_alpha(0.7)
                line_artist.set_linewidth(1.5)

        if self.highlight_annotation:
            self.highlight_annotation.remove()
            self.highlight_annotation = None

        self.selected_exp_id = None
        self.canvas.draw()
    
    def update_plot_mode3(self):
        """Update plot for Mode 3 (Photon Type)"""
        photon_type = self.photon_type_combo.get()
        cell_filter = self.photon_cell_filter_combo.get()
        
        if not photon_type:
            return
        
        # Get experiments
        experiments = self.data_manager.get_experiments_by_photon(
            photon_type, cell_filter if cell_filter != 'All' else None)
        
        # Filter to unique photon references only
        # Keep first experiment for each (PublicationID, PhotonExpNum) combination
        # that has actual raw data
        seen_refs = set()
        unique_experiments = []
        for _, exp in experiments.iterrows():
            pub_id = exp['PublicationID']
            photon_exp = int(exp['PhotonExpNum']) if pd.notna(exp['PhotonExpNum']) else 0
            ref_key = (pub_id, photon_exp)
            # Only include if we haven't seen this ref AND it has raw data
            if ref_key not in seen_refs:
                if self.data_manager.get_photon_raw_data(pub_id, photon_exp) is not None:
                    seen_refs.add(ref_key)
                    unique_experiments.append(exp)
        
        self.current_experiments = pd.DataFrame(unique_experiments)
        
        # Update statistics
        n_exp = len(self.current_experiments)
        cells_in_results = self.current_experiments['Cells'].unique() if not self.current_experiments.empty else []
        
        self.stats_label.config(text=
            f'Photon: {photon_type}\n'
            f'Unique curves: {n_exp}\n'
            f'Cell lines: {len(cells_in_results)}'
        )
        
        # Update table
        self.update_table()
        
        # Update plot
        self.plot_survival_data_mode3()
        
        self.status_var.set(f'Displaying {n_exp} unique photon curves for {photon_type}')
    
    def plot_survival_data_mode3(self):
        """Plot survival curves for Mode 3 (Photon radiation)
        
        Note: In photon mode, multiple experiments may share the same photon reference data
        (same PublicationID + PhotonExpNum). Each unique photon reference is plotted once,
        and all experiments sharing that reference are linked to the same plot element.
        """
        self.ax.clear()
        self.plot_elements = {}
        self.selected_exp_id = None
        self.highlight_annotation = None
        
        experiments = self.current_experiments
        
        if experiments.empty:
            self.setup_plot()
            self.ax.text(0.5, 0.5, 'No data to display', 
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize=14, color='gray')
            self.canvas.draw()
            return
        
        # Color map for different cell lines
        cell_lines = experiments['Cells'].dropna().unique()
        colors = plt.cm.tab10(np.linspace(0, 1, max(1, len(cell_lines))))
        cell_color_map = dict(zip(sorted(cell_lines), colors))
        
        show_raw = self.show_raw_var.get()
        show_lq = self.show_lq_var.get()
        
        plotted_cells = set()
        plotted_photon_refs = {}  # {(pub_id, photon_exp): {'scatter': artist, 'line': artist}}
        
        for idx, (_, exp) in enumerate(experiments.iterrows()):
            exp_id = exp['ExpID']
            cell = exp['Cells']
            pub_id = exp['PublicationID']
            photon_exp = int(exp['PhotonExpNum']) if pd.notna(exp['PhotonExpNum']) else 0
            photon_ref_key = (pub_id, photon_exp)
            color = cell_color_map.get(cell, 'blue')
            
            # Check if this photon reference was already plotted
            if photon_ref_key in plotted_photon_refs:
                # Link this experiment to existing plot elements
                self.plot_elements[exp_id] = plotted_photon_refs[photon_ref_key].copy()
                continue
            
            # Get photon raw data
            raw_data = self.data_manager.get_photon_raw_data(pub_id, photon_exp)
            has_raw = raw_data is not None
            
            # Create label only for first occurrence of each cell line
            label_raw = f'{cell}' if cell not in plotted_cells else None
            label_lq = None
            
            # Initialize storage for this photon reference
            plotted_photon_refs[photon_ref_key] = {}
            
            # Plot raw data points
            if show_raw and has_raw:
                scatter_artist = self.ax.scatter(
                    raw_data['doses'], raw_data['sf'],
                    color=color, s=50, alpha=0.8, zorder=3,
                    label=label_raw
                )
                plotted_photon_refs[photon_ref_key]['scatter'] = scatter_artist
                self.plot_elements.setdefault(exp_id, {})['scatter'] = scatter_artist
                plotted_cells.add(cell)
            
            # Plot LQ curve using photon parameters (ax_paper or ax_fit)
            # Try ax_paper first, then fall back to ax_fit
            alpha_val = exp['ax_paper'] if pd.notna(exp['ax_paper']) else exp.get('ax_fit')
            beta_col = 'bx_paper' if pd.notna(exp['ax_paper']) else 'bx_fit'
            
            if show_lq and pd.notna(alpha_val):
                alpha = alpha_val
                beta = exp[beta_col] if pd.notna(exp.get(beta_col)) else 0
                
                doses = np.linspace(0, 12, 100)
                sf = np.exp(-(alpha * doses + beta * doses**2))
                
                linestyle = '-' if has_raw else '--'
                linewidth = 1.5 if has_raw else 2
                
                # Label for LQ curve if no raw data was plotted for this cell
                if cell not in plotted_cells:
                    label_lq = f'{cell}'
                    plotted_cells.add(cell)
                
                line_artist, = self.ax.plot(
                    doses, sf, color=color, linestyle=linestyle,
                    linewidth=linewidth, alpha=0.7, zorder=2, label=label_lq
                )
                plotted_photon_refs[photon_ref_key]['line'] = line_artist
                self.plot_elements.setdefault(exp_id, {})['line'] = line_artist

        
        # Set up axes
        self.ax.set_xlabel('Dose (Gy)', fontsize=12)
        self.ax.set_ylabel('Survival Fraction', fontsize=12)
        self.ax.set_yscale('log')
        self.ax.set_xlim(0, 12)
        self.ax.set_ylim(1e-4, 1.5)
        self.ax.grid(True, alpha=0.3)
        
        # Title
        photon_type = self.photon_type_combo.get()
        self.ax.set_title(f'Cell Survival with {photon_type} Photon Radiation', fontsize=14)
        
        # Legend
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            if len(handles) > 12:
                handles = handles[:12]
                labels = labels[:12]
            self.ax.legend(handles, labels, loc='lower left', fontsize=8, ncol=2)
        
        self.fig.tight_layout()
        self.canvas.draw()
    
    def update_plot_mode2(self):
        """Update plot for Mode 2 (Radiation Type → LET)"""
        rad_type = self.rad_type_combo.get()
        let_str = self.let_combo.get()
        cell_filter = self.cell_filter_combo.get()
        
        if not rad_type:
            return
        
        # Parse LET value
        let_val = None
        if let_str and let_str != 'All':
            let_val = float(let_str)
        
        # Get experiments
        self.current_experiments = self.data_manager.get_experiments_by_ion_let(
            rad_type, let_val, cell_filter if cell_filter != 'All' else None)
        
        # Update statistics
        n_exp = len(self.current_experiments)
        n_raw = sum(1 for _, exp in self.current_experiments.iterrows() 
                   if self.data_manager.get_ion_raw_data(exp['ExpID']) is not None)
        
        cells_in_results = self.current_experiments['Cells'].unique() if not self.current_experiments.empty else []
        
        self.stats_label.config(text=
            f'Ion: {rad_type}\n'
            f'LET: {let_str} keV/μm\n'
            f'Experiments: {n_exp}\n'
            f'With raw data: {n_raw}\n'
            f'Cell lines: {len(cells_in_results)}'
        )
        
        # Update table
        self.update_table()
        
        # Update plot
        self.plot_survival_data_mode2()
        
        self.status_var.set(f'Displaying {n_exp} experiments for {rad_type} (LET={let_str})')
    
    def plot_survival_data_mode2(self):
        """Plot survival curves for Mode 2 (grouped by cell line instead of LET)"""
        self.ax.clear()
        self.plot_elements = {}
        self.selected_exp_id = None
        self.highlight_annotation = None
        
        experiments = self.current_experiments
        
        if experiments.empty:
            self.setup_plot()
            self.ax.text(0.5, 0.5, 'No data to display', 
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize=14, color='gray')
            self.canvas.draw()
            return
        
        # Color map for different cell lines
        cell_lines = experiments['Cells'].dropna().unique()
        colors = plt.cm.tab10(np.linspace(0, 1, max(1, len(cell_lines))))
        cell_color_map = dict(zip(sorted(cell_lines), colors))
        
        show_raw = self.show_raw_var.get()
        show_lq = self.show_lq_var.get()
        
        plotted_cells = set()
        
        for idx, (_, exp) in enumerate(experiments.iterrows()):
            exp_id = exp['ExpID']
            cell = exp['Cells']
            let = exp['LET'] if pd.notna(exp['LET']) else 0
            color = cell_color_map.get(cell, 'blue')
            
            # Get raw data
            raw_data = self.data_manager.get_ion_raw_data(exp_id)
            has_raw = raw_data is not None
            
            # Create label only for first occurrence of each cell line
            label_raw = f'{cell} (LET={let:.1f})' if cell not in plotted_cells else None
            label_lq = None
            
            # Plot raw data points
            if show_raw and has_raw:
                scatter_artist = self.ax.scatter(
                    raw_data['doses'], raw_data['sf'],
                    color=color, s=50, alpha=0.8, zorder=3,
                    label=label_raw
                )
                self.plot_elements.setdefault(exp_id, {})['scatter'] = scatter_artist
                plotted_cells.add(cell)
            
            # Plot LQ curve
            if show_lq and pd.notna(exp['ai_paper']):
                alpha = exp['ai_paper']
                beta = exp['bi_paper'] if pd.notna(exp['bi_paper']) else 0
                
                doses = np.linspace(0, 12, 100)
                sf = np.exp(-(alpha * doses + beta * doses**2))
                
                linestyle = '-' if has_raw else '--'
                linewidth = 1.5 if has_raw else 2
                
                # Label for LQ curve if no raw data was plotted for this cell
                if cell not in plotted_cells:
                    label_lq = f'{cell} (LET={let:.1f})'
                    plotted_cells.add(cell)
                
                line_artist, = self.ax.plot(
                    doses, sf, color=color, linestyle=linestyle,
                    linewidth=linewidth, alpha=0.7, zorder=2, label=label_lq
                )
                self.plot_elements.setdefault(exp_id, {})['line'] = line_artist
        
        # Set up axes
        self.ax.set_xlabel('Dose (Gy)', fontsize=12)
        self.ax.set_ylabel('Survival Fraction', fontsize=12)
        self.ax.set_yscale('log')
        self.ax.set_xlim(0, 12)
        self.ax.set_ylim(1e-4, 1.5)
        self.ax.grid(True, alpha=0.3)
        
        # Title
        rad_type = self.rad_type_combo.get()
        let_str = self.let_combo.get()
        self.ax.set_title(f'{rad_type} Ion Survival Curves (LET={let_str} keV/μm)', fontsize=14)
        
        # Legend
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            if len(handles) > 12:
                handles = handles[:12]
                labels = labels[:12]
            self.ax.legend(handles, labels, loc='lower left', fontsize=8, ncol=2)
        
        self.fig.tight_layout()
        self.canvas.draw()
    
    def update_plot(self):
        """Update the survival curve plot"""
        # Check mode
        mode = self.selection_mode.get()
        if mode == 'radiation_let':
            self.update_plot_mode2()
            return
        elif mode == 'photon':
            self.update_plot_mode3()
            return
        
        cell_line = self.cell_combo.get()
        ion = self.ion_combo.get()
        
        if not cell_line or not ion:
            return
        
        # Get experiments
        self.current_experiments = self.data_manager.get_experiments(cell_line, ion)
        
        # Update statistics
        n_exp = len(self.current_experiments)
        n_raw = sum(1 for _, exp in self.current_experiments.iterrows() 
                   if self.data_manager.get_ion_raw_data(exp['ExpID']) is not None)
        
        let_range = ''
        if not self.current_experiments.empty and 'LET' in self.current_experiments.columns:
            let_vals = self.current_experiments['LET'].dropna()
            if len(let_vals) > 0:
                let_range = f'\nLET range: {let_vals.min():.1f} - {let_vals.max():.1f} keV/μm'
        
        self.stats_label.config(text=
            f'Cell: {cell_line}\n'
            f'Ion: {ion}\n'
            f'Experiments: {n_exp}\n'
            f'With raw data: {n_raw}'
            f'{let_range}'
        )
        
        # Update table
        self.update_table()
        
        # Update plot
        self.plot_survival_data()
        
        self.status_var.set(f'Displaying {n_exp} experiments for {cell_line} + {ion}')
    
    def plot_survival_data(self):
        """Plot survival curves for current experiments"""
        self.ax.clear()
        self.plot_elements = {}
        self.selected_exp_id = None
        self.highlight_annotation = None
        
        experiments = self.current_experiments
        
        if experiments.empty:
            self.setup_plot()
            self.ax.text(0.5, 0.5, 'No data to display', 
                        transform=self.ax.transAxes, ha='center', va='center',
                        fontsize=14, color='gray')
            self.canvas.draw()
            return
        
        # Color map for different LET values
        let_values = experiments['LET'].dropna().unique()
        colors = plt.cm.viridis(np.linspace(0, 0.9, max(1, len(let_values))))
        let_color_map = dict(zip(sorted(let_values), colors))
        
        show_raw = self.show_raw_var.get()
        show_lq = self.show_lq_var.get()
        show_photon = self.show_photon_var.get()
        
        for idx, (_, exp) in enumerate(experiments.iterrows()):
            exp_id = exp['ExpID']
            let = exp['LET'] if pd.notna(exp['LET']) else 0
            color = let_color_map.get(let, 'blue')
            
            # Get raw data
            raw_data = self.data_manager.get_ion_raw_data(exp_id)
            has_raw = raw_data is not None
            
            # Plot raw data points
            if show_raw and has_raw:
                scatter_artist = self.ax.scatter(
                    raw_data['doses'], raw_data['sf'],
                    color=color, s=50, alpha=0.8, zorder=3,
                    label=f'ExpID {exp_id} (LET={let:.1f})'
                )
                self.plot_elements.setdefault(exp_id, {})['scatter'] = scatter_artist
            
            # Plot LQ curve
            if show_lq and pd.notna(exp['ai_paper']):
                alpha = exp['ai_paper']
                beta = exp['bi_paper'] if pd.notna(exp['bi_paper']) else 0
                
                doses = np.linspace(0, 12, 100)
                sf = np.exp(-(alpha * doses + beta * doses**2))
                
                linestyle = '-' if has_raw else '--'
                linewidth = 1.5 if has_raw else 2
                
                label = f'LQ (α={alpha:.3f}, β={beta:.4f})' if not has_raw else None
                line_artist, = self.ax.plot(
                    doses, sf, color=color, linestyle=linestyle,
                    linewidth=linewidth, alpha=0.7, zorder=2, label=label
                )
                self.plot_elements.setdefault(exp_id, {})['line'] = line_artist
        
        # Plot photon reference if requested
        if show_photon and not experiments.empty:
            first_exp = experiments.iloc[0]
            pub_id = first_exp['PublicationID']
            photon_exp = first_exp['PhotonExpNum']
            
            photon_raw = self.data_manager.get_photon_raw_data(pub_id, photon_exp)
            if photon_raw is not None:
                self.ax.scatter(photon_raw['doses'], photon_raw['sf'], 
                               color='black', marker='x', s=40, alpha=0.6, zorder=2,
                               label='Photon ref (raw)')
            
            # Plot photon LQ curve
            if pd.notna(first_exp['ax_paper']):
                alpha_x = first_exp['ax_paper']
                beta_x = first_exp['bx_paper'] if pd.notna(first_exp['bx_paper']) else 0
                doses = np.linspace(0, 12, 100)
                sf_x = np.exp(-(alpha_x * doses + beta_x * doses**2))
                self.ax.plot(doses, sf_x, 'k--', linewidth=1.5, alpha=0.5, 
                            label='Photon ref (LQ)')
        
        # Set up axes
        self.ax.set_xlabel('Dose (Gy)', fontsize=12)
        self.ax.set_ylabel('Survival Fraction', fontsize=12)
        self.ax.set_yscale('log')
        self.ax.set_xlim(0, 12)
        self.ax.set_ylim(1e-4, 1.5)
        self.ax.grid(True, alpha=0.3)
        
        # Title
        if not experiments.empty:
            cell = experiments.iloc[0]['Cells']
            ion = experiments.iloc[0]['Ion']
            self.ax.set_title(f'{cell} Cell Survival with {ion} Ions', fontsize=14)
        
        # Legend
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:
            # Limit legend entries
            if len(handles) > 10:
                handles = handles[:10]
                labels = labels[:10]
            self.ax.legend(handles, labels, loc='lower left', fontsize=8, ncol=2)
        
        self.fig.tight_layout()
        self.canvas.draw()
    
    def update_table(self):
        """Update the experiment details table"""
        # Clear existing items
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        for _, exp in self.current_experiments.iterrows():
            # ExpID
            exp_id = int(exp['ExpID'])
            
            # Cell line
            cell = str(exp['Cells']) if pd.notna(exp['Cells']) else 'N/A'
            
            # Publication
            pub = str(exp['PublicationName'])
            
            # LET
            let_str = f"{exp['LET']:.2f}" if pd.notna(exp['LET']) else 'N/A'
            
            # Energy
            energy_str = f"{exp['Energy']:.1f}" if pd.notna(exp['Energy']) else 'N/A'
            
            # Alpha
            alpha_str = f"{exp['ai_paper']:.4f}" if pd.notna(exp['ai_paper']) else 'N/A'
            
            # Beta
            beta_str = f"{exp['bi_paper']:.5f}" if pd.notna(exp['bi_paper']) else 'N/A'
            
            # D10
            if pd.notna(exp['ai_paper']):
                alpha = exp['ai_paper']
                beta = exp['bi_paper'] if pd.notna(exp['bi_paper']) else 0
                target = np.log(10)
                if beta == 0:
                    d10 = target / alpha if alpha > 0 else 0
                else:
                    discriminant = alpha**2 + 4*beta*target
                    d10 = (-alpha + np.sqrt(discriminant)) / (2*beta) if discriminant > 0 else 0
                d10_str = f"{d10:.2f}"
            else:
                d10_str = 'N/A'
            
            # Photon reference - shows which experiments share the same photon reference data
            mode = self.selection_mode.get()
            if mode == 'photon':
                pub_id = exp['PublicationID']
                photon_exp = int(exp['PhotonExpNum']) if pd.notna(exp['PhotonExpNum']) else 0
                photon_ref = f"P{pub_id}-{photon_exp}"
                has_raw = 'Yes' if self.data_manager.get_photon_raw_data(pub_id, photon_exp) is not None else 'No'
            else:
                photon_ref = '-'
                has_raw = 'Yes' if self.data_manager.get_ion_raw_data(exp['ExpID']) is not None else 'No'
            
            self.tree.insert('', tk.END, values=(exp_id, cell, pub, let_str, energy_str, 
                                                  alpha_str, beta_str, d10_str, photon_ref, has_raw))
    
    def clear_plot(self):
        """Clear the plot"""
        self.ax.clear()
        self.setup_plot()
        self.canvas.draw()
        
        # Clear table
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        self.stats_label.config(text='No data displayed')
    
    def save_plot(self):
        """Save the current plot"""
        from tkinter import filedialog
        
        filename = filedialog.asksaveasfilename(
            defaultextension='.png',
            filetypes=[('PNG files', '*.png'), ('PDF files', '*.pdf'), ('All files', '*.*')],
            title='Save Plot'
        )
        
        if filename:
            self.fig.savefig(filename, dpi=150, bbox_inches='tight')
            self.status_var.set(f'Plot saved to {filename}')


def main():
    # Determine PIDE directory
    if len(sys.argv) > 1:
        pide_dir = sys.argv[1]
    else:
        # Default paths to try
        script_dir = Path(__file__).parent
        possible_paths = [
            script_dir / 'PIDE3.4',
            script_dir.parent / 'PIDE3.4',
            Path('/home/user/Geant4Projects/CellStateTransitionPaper/RADCellSimulation/PIDE3.4'),
        ]
        
        pide_dir = None
        for path in possible_paths:
            if path.exists():
                pide_dir = str(path)
                break
        
        if pide_dir is None:
            print("Error: Could not find PIDE3.4 folder")
            print("Usage: python3 pide_gui.py [path_to_PIDE3.4_folder]")
            sys.exit(1)
    
    print(f"Using PIDE data from: {pide_dir}")
    
    # Create application
    root = tk.Tk()
    
    # Set theme
    style = ttk.Style()
    style.theme_use('clam')  # Use 'clam' theme for better appearance
    
    # Create GUI
    print("Creating GUI...")
    app = PIDEExplorerGUI(root, pide_dir)
    
    # Force window to front
    root.lift()
    root.attributes('-topmost', True)
    root.after(100, lambda: root.attributes('-topmost', False))
    root.focus_force()
    
    print("GUI created. Starting main loop...")
    print("If you don't see the window, check other workspaces/monitors.")
    
    # Run
    root.mainloop()


if __name__ == '__main__':
    main()
