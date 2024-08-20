#!/usr/bin/env python3

################################################################
# Program: hDNApipe-GUI
# Version 1.0
# Author: Zhang Yaxin (zhangyaxin@mail.ustc.edu.cn)
#################################################################

import os
import tkinter as tk
import tkinter.ttk as ttk
import tkinter.font as tkFont
import subprocess
import threading
import time
from tkinter import filedialog
from tkinter import messagebox

## path
dnapipe_dir = os.path.dirname(os.path.abspath(__file__))

## style

# common
style_button_big = { "bg": "skyblue", "foreground": "white", "font": ("Nimbus Roman", 15, "bold"), "width": 12}
style_button_browse = { "bg": "steelblue", "foreground": "white", "width": 6}
style_entry = {"bg": "whitesmoke"}
style_label_header = {"bg": "white", "anchor": "w", "font": ("Default", 15)}
# adv
style_group= {"bg": "white", "padx": 10, "pady": 5, "labelanchor": "nw"}
style_radiobutton = { "bg": "white", "highlightthickness":0 , "anchor": "w", "width": 15 }
style_checkbox = { "bg": "white", "highlightthickness":0 , "anchor": "w", "width": 20 }
style_label = {"bg": "white", "anchor": "w"}
# plot
pos_label_header = {"sticky":'w', "padx": 30, "pady": 10, "columnspan": 2}
pos_label_col1st = {"sticky":'w', "padx": (60,10), "pady": 2}

## browse functions

def select_file(entry):
    filename=filedialog.askopenfilename()
    entry.delete(0,tk.END)
    entry.insert(0,filename)
    entry.xview_moveto(1)

def select_folder(entry):
    foldername = filedialog.askdirectory()
    entry.delete(0,tk.END)
    entry.insert(0,foldername)
    entry.xview_moveto(1)

## GUI

class init_page(tk.Frame):

    def __init__(self, parent):
        super().__init__(parent)
        self.configure(bg="white")
        image_path = dnapipe_dir + r"/images/logo_basic.png"
        self.imageobj = tk.PhotoImage(file=image_path)
        label_img = tk.Label(self, image = self.imageobj, anchor='center', bg="white")
        label_img.pack(side='top', anchor='center', expand=True, pady=(10,0))
        
        group = tk.LabelFrame(self, text="",borderwidth=0, **style_group)

        text1 = "For first-time use initiation:"
        text2 = "Step 1. Click \"Initiate\" button to install necessary resources."
        text3 = "Step 2. Click \"Reference\" button to set or get hg38 reference genome."
        text4 = "* For very detailed information, please visit https://github.com/TJ-lab-ustc/hDNApipe"
        for i in range(1,5):
            text="text"+str(i)
            message_text = locals()[text]
            tk.Message(group, text=message_text, anchor='center', width=500, bg="white", font=("TkDefaultFont", 15)).grid(row=i, column=0, columnspan=5, pady=(0, 5), sticky='w')

        tk.Button(group, text="Initiate", command=self.run_init, **style_button_big).grid(row=99, column=1, sticky='ew', pady=(15,0))
        tk.Button(group, text="Reference", command=self.set_ref, **style_button_big).grid(row=99, column=3, sticky='ew', pady=(15,0))
        
        group.pack(side='bottom', anchor='center', expand=True)
        
    def run_init(self):
        
        def run_init_shell():
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            while True:
                output = proc.stdout.readline().decode().strip()
                if not output and proc.poll() is not None:
                    break
                output_box.config(state="normal")
                output_box.insert(tk.END, output + "\n")
                output_box.config(state="disabled")
                output_box.see(tk.END)

            proc.communicate()
            if (proc.returncode != 0):
                tk.messagebox.showinfo("INFO", "Initiation failed. See log in the output window.")

        output_box.config(state="normal")
        output_box.insert(tk.END, "Initiating..." + "\n")
        output_box.config(state="disabled")
        dnapipe_script = dnapipe_dir + "/dnapipe"
        command = ["bash", dnapipe_script, "init", "--annot"]
        t = threading.Thread(target=run_init_shell)
        t.start()

    def set_ref(self):

        def run_set_ref():

            def set_ref_shell():
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while True:
                    output = proc.stdout.readline().decode().strip()
                    if not output and proc.poll() is not None:
                        break
                    output_box.config(state="normal")
                    output_box.insert(tk.END, output + "\n")
                    output_box.config(state="disabled")
                    output_box.see(tk.END)
                proc.communicate()
                if (proc.returncode != 0):
                    tk.messagebox.showinfo("INFO", "Execution failed. See log in output window.")

            if ref_option.get() == "download":
                path = dnapipe_dir
            else:
                path = ref_input.get()
            top.destroy()
            output_box.config(state="normal")
            output_box.insert(tk.END, "Reference setting..." + "\n")
            output_box.config(state="disabled")
            para = "--" + ref_option.get()
            dnapipe_script = dnapipe_dir + "/dnapipe"
            command = ["bash", dnapipe_script, "ref", para, path]
            t = threading.Thread(target=set_ref_shell)
            t.start()

        top = tk.Toplevel(bg="white")
        top.geometry("300x180")
        top.title("Set reference")
        top.wm_transient(self)
        tk.Label(top, text="How to set your reference?", bg="white").grid(row=0, column=0, sticky='w', columnspan=2)
        ref_option = tk.StringVar()
        ref_option.set("download")
        select_path = tk.StringVar()
        style_radiobutton1 = { "bg": "white", "highlightthickness":0 , "anchor": "w", "width": 21 }
        tk.Radiobutton(top, text="Automatic download", variable=ref_option, value="download", **style_radiobutton1).grid(row=1 ,column=0, sticky='w', columnspan=2)
        tk.Radiobutton(top, text="Choose from files", variable=ref_option, value="set", **style_radiobutton1).grid(row=2, column=0, sticky='w', columnspan=2)
        ref_input = tk.Entry(top, width=25, textvariable=select_path)
        ref_input.grid(row=3, column=0, pady=10, padx=(15,0))
        tk.Button(top, text="Browse", command=lambda: select_file(ref_input), **style_button_browse).grid(row=3, column=1, pady=(10,0), padx=(5,0))
        tk.Button(top, text="OK", command=run_set_ref, **style_button_big).grid(row=99 ,column=0, pady=10)

class advance_page(tk.Frame):

    def __init__(self, parent):
        super().__init__(parent)
        self.configure(bg="white")
        self.var_group = tk.LabelFrame(self, text="",borderwidth=0, **style_group)
        self.align_pre_wigets()
        self.var_wigets()
        self.annot_wigets()
        self.computing_wigets()
        self.set_adv_option_default()
        self.get_advance()
        tk.Button(self.var_group, text="Save", command=self.adv_run_button_click, **style_button_big).grid(column=0, row=4, pady=10)
        tk.Button(self.var_group, text="Default", command=self.adv_default_button_click, **style_button_big).grid(column=1, row=4, pady=10)
        self.var_group.pack(anchor='center', expand=True)

    def var_wigets(self):
        group = tk.LabelFrame(self.var_group, text="Variants calling", **style_group)
        label_list = ["SNV/INDEL calling tools:", "SV calling tools:", "CNV calling tools:", "Else:"]
        for i in range(len(label_list)):
            tk.Label(group, text=label_list[i], **style_label).grid(row=4*i, column=0, sticky='w', columnspan=2, pady=1)
        var_tool_list = ["Strelka", "DeepVariant", "GATK", "Manta", "Lumpy", "DELLY", "CNVkit","DELLY"]
        self.vartool_checkbutton = {}
        self.vartool_variable = {}
        for i in range(len(var_tool_list)):
            self.vartool_variable[i] = tk.BooleanVar()
            self.vartool_checkbutton[i] = tk.Checkbutton(group, text=var_tool_list[i], variable=self.vartool_variable[i], **style_checkbox)
            self.vartool_checkbutton[i].grid(row=i+i//3+1, column=0, sticky='w', padx=10, columnspan=2, pady=2)

        tk.Label(group, text=" -- CNV bin size(bp):", **style_label).grid(row=14, column=0, sticky='w', padx=(10,5), pady=1)
        self.cnv_bin_entry = tk.Entry(group, width=8, **style_entry)
        self.cnv_bin_entry.grid(row=14, column=1, sticky='w', pady=1)
        group.grid(row=2, column=0, rowspan=2, sticky='ew', padx=(0, 10))

    def align_pre_wigets(self):
        group = tk.LabelFrame(self.var_group, text="Align and preprocess", **style_group)
        pre_lists = {0: "Force Replacing RG" ,1: "Remove Dupicates", 2: "Soft Clipping"}
        self.pre_variable = {}
        self.pre_checkbutton = {}
        for i in range(len(pre_lists)):
            self.pre_variable[i] = tk.BooleanVar()
            self.pre_checkbutton[i] = tk.Checkbutton(group, text=pre_lists[i], variable=self.pre_variable[i], **style_checkbox)
            self.pre_checkbutton[i].grid(row=i, column=0, padx=(0,10), sticky='w', columnspan=2, pady=1)
        tk.Label(group, text=" -- MAPQ threshold:", **style_label, width=17).grid(row=3, column=0, sticky='w', pady=1)
        self.mapq_entry = tk.Entry(group, width=5, **style_entry)
        self.mapq_entry.grid(row=3, column=1, sticky='w', pady=1)
        group.grid(row=2, column=1, columnspan=2, sticky='ew')

    def annot_wigets(self):

        def click_annot_all():
            for i in range(len(list)):
                self.annot_variable[i].set("True")
        
        def clear_annot_all():
            for i in range(len(list)):
                self.annot_variable[i].set("False")

        group = tk.LabelFrame(self.var_group, text="Annotation", **style_group)
        list = ["Default annotation", "Gene symbol", "Check existing", "Pathogenicity", "Group frequency", "Clinvar", "Report"]
        self.annot_variable = {}
        self.annot_checkbutton = {}
        for i in range(len(list)):
            self.annot_variable[i] = tk.BooleanVar()
            self.annot_checkbutton[i] = tk.Checkbutton(group, text=list[i], variable=self.annot_variable[i], **style_checkbox)
            self.annot_checkbutton[i].grid(row=i+6, column=0, padx=(0,10), sticky='w', columnspan=2, pady=1)
        tk.Button(group, text = "Click ALL", command=click_annot_all, bg="whitesmoke").grid(row=20, column=0)
        tk.Button(group, text = "Clear", command=clear_annot_all, bg="whitesmoke").grid(row=20, column=1)
        group.grid(row=3, column=1, columnspan=2, sticky='ew')
        
    def computing_wigets(self):
        group = tk.LabelFrame(self.var_group, text="Computing paramters", **style_group)
        tk.Label(group, text="Threads:", **style_label).grid(row=1, column=0, padx=(0,10), sticky='w')
        self.thread_entry = tk.Entry(group, width=10, **style_entry)
        self.thread_entry.grid(row=1, column=1, sticky='w')
        tk.Label(group, text="Memory(G):", **style_label).grid(row=1, column=2, padx=(110,10), sticky='w')
        self.memory_entry = tk.Entry(group, width=10, **style_entry)
        self.memory_entry.grid(row=1, column=3, sticky='w', padx=(0, 30))
        group.grid(row=1, column=0, columnspan=2, sticky='ew', pady=(0, 10))

    def set_adv_option_default(self):
        for i in range(len(self.vartool_variable)):
            if i in [0, 3, 6]:
                self.vartool_variable[i].set("True")
                True
            else:
                self.vartool_variable[i].set("False")
        for i in range(len(self.pre_variable)):
            self.pre_variable[i].set("False")
        for i in range(len(self.annot_variable)):
            self.annot_variable[i].set("False")
        entry_list = [self.thread_entry, self.memory_entry, self.mapq_entry, self.cnv_bin_entry]
        for entry in entry_list:
            entry.delete(0,tk.END)
        total_thread = os.cpu_count()
        thread_default = round(total_thread/2) if round(total_thread/2) >= 1 else 1
        self.thread_entry.insert(0, thread_default)
        total_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        memory_default = round(total_memory/1024/1024/1024/2) if round(total_memory/2) >= 1 else 1
        self.memory_entry.insert(0, memory_default)
        self.mapq_entry.insert(0, "30")
        self.pre_variable[1].set("True")
        self.cnv_bin_entry.insert(0, "10000")
    
    def get_advance(self):
        short_list = {
            "strelka": self.vartool_variable[0].get(),
            "deepvariant": self.vartool_variable[1].get(),
            "gatk": self.vartool_variable[2].get()
        }
        sv_list = {
            "manta": self.vartool_variable[3].get(),
            "lumpy": self.vartool_variable[4].get(),
            "delly": self.vartool_variable[5].get()
        }
        cnv_list = {
            "cnvkit": self.vartool_variable[6].get(),
            "delly_cnv": self.vartool_variable[7].get(),
        }
        short_tools = ','.join([key for key, value in short_list.items() if value])
        sv_tools = ','.join([key for key, value in sv_list.items() if value])
        cnv_tools = ','.join([key for key, value in cnv_list.items() if value])

        global advance_option
        advance_option = {
            "threads": self.thread_entry.get(),
            "memory": self.memory_entry.get(),
            
            "short_tools": short_tools,
            "sv_tools": sv_tools,
            "cnv_tools": cnv_tools,
            "bin": self.cnv_bin_entry.get(),

            "forceRG": self.pre_variable[0].get(),
            "rmdup": self.pre_variable[1].get(),
            "softclip": self.pre_variable[2].get(),
            "mapq": self.mapq_entry.get(),

            # annotation
            "default": self.annot_variable[0].get(),
            "symbol": self.annot_variable[1].get(),
            "exist": self.annot_variable[2].get(),
            "pathogenicity": self.annot_variable[3].get(),
            "frequency": self.annot_variable[4].get(),
            "clinvar": self.annot_variable[5].get(),
            "report": self.annot_variable[6].get()
        }
        if advance_option["threads"] == "":
            total_thread = os.cpu_count()
            advance_option["threads"] = round(total_thread/2) if round(total_thread/2) >= 1 else 1
        advance_option["threads"] = str(advance_option["threads"])
        if advance_option["memory"] == "":
            total_memory = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
            advance_option["memory"] = round(total_memory/1024/1024/1024/2) if round(total_memory/2) >= 1 else 1
        advance_option["memory"] = str(advance_option["memory"]) + "G"
        if not advance_option["mapq"]:
            advance_option["mapq"] = "30"
        if not short_tools:
            advance_option["short_tools"] = "strelka"
        if not sv_tools:
            advance_option["sv_tools"] = "manta"
        if not cnv_tools:
            advance_option["cnv_tools"] = "cnvkit"
        
    def print_adv_change(self):
        output_box.config(state="normal")
        output_box.insert(tk.END, "Advanced options has changed." + "\n" + str(advance_option) + "\n" + "\n")
        output_box.config(state="disabled")
        output_box.see(tk.END)

    def adv_run_button_click(self): 
        return [self.get_advance(), self.print_adv_change()]
    
    def adv_default_button_click(self):
        return [self.set_adv_option_default(), self.get_advance(), self.print_adv_change()]

class var_page(tk.Frame):

    def __init__(self, parent):
        super().__init__(parent)

        self.configure(bg="white")
        self.create_option_widgets()

    def create_option_widgets(self):

        group = tk.LabelFrame(self, text="",borderwidth=0, **style_group)

        style_radiobutton = { "bg": "white", "highlightthickness":0 , "anchor": "w", "width": 5 }
        style_checkbox = { "bg": "white", "highlightthickness":0 , "anchor": "w", "width": 5 }
        grid_label_header = {"sticky": "w", "pady": 10, "columnspan": 2}
        grid_label_col1st = {"sticky": "w", "padx": (0, 10), "pady": 7}

        tk.Label(group, text="Analysis Options", **style_label_header).grid(row=0, column=0, sticky='w', columnspan=2)

        row1 = 1
        label_list2 = ["Detection Mode:", "Sequencing Method:", "Input file Type:", "Variant calling:", "", "Region(.bed)*:", "Output directory:"]
        for i in range(len(label_list2)):
            tk.Label(group, text=label_list2[i], **style_label).grid(row=row1+i, column=0, **grid_label_col1st)
        
        detect_mode_lists = ["germline", "somatic"]
        self.detect_mode_combobox = ttk.Combobox(group, values=detect_mode_lists, justify="center", state="readonly", takefocus=0)
        self.detect_mode_combobox.grid(row=1, column=1, columnspan=2, sticky='w')
        row1+=1

        self.seq_method = tk.StringVar()
        tk.Radiobutton(group, text="WGS", variable=self.seq_method, value="wgs", **style_radiobutton).grid(row=row1, column=1, sticky='w')
        wes_radiobutton = tk.Radiobutton(group, text="WES/Target", variable=self.seq_method, value="wes", **style_radiobutton)
        wes_radiobutton.config(width=10)
        wes_radiobutton.grid(row=row1, column=2, sticky='w', columnspan=2)
        row1+=1
        
        self.file_type = tk.StringVar()
        type_list = {0:"fastq", 1:"bam"}
        for i in range(len(type_list)):
            type=type_list[i]
            tk.Radiobutton(group, text=type, variable=self.file_type, value=type, **style_radiobutton).grid(row=row1, column=i+1, sticky='w')
        row1+=1

        var_lists = {0: "SNV & INDEL", 1: "SV", 2: "CNV"}
        self.vartype_variable = {}
        self.vartype_checkbutton = {}
        for i in range(len(var_lists)):
            self.vartype_variable[i] = tk.BooleanVar()
            self.vartype_checkbutton[i] = tk.Checkbutton(group, text=var_lists[i], variable=self.vartype_variable[i], **style_checkbox)
            if i==0:
                self.vartype_checkbutton[i].grid(row=row1, column=i + 1, sticky='w', columnspan=2)
                self.vartype_checkbutton[i].config(width=12)
            else:
                self.vartype_checkbutton[i].grid(row=row1+1, column=i, sticky='w')
        row1+=2
        
        self.entry_region = tk.Entry(group, width=30, **style_entry)
        self.entry_region.grid(row=row1, column=1, columnspan=2, sticky='w')
        tk.Button(group, text="Browse", command=lambda: select_file(self.entry_region), **style_button_browse).grid(row=row1, column=3, sticky='w', padx=(10,0))
        row1+=1

        self.entry_outdir = tk.Entry(group, width=30, **style_entry)
        self.entry_outdir.grid(row=row1, column=1, columnspan=2, sticky='w')
        tk.Button(group, text="Browse", command=lambda: select_folder(self.entry_outdir), **style_button_browse).grid(row=row1, column=3, sticky='w', padx=(10,0))
        row1+=1

        row1+=len(label_list2)
        tk.Label(group, text="Sample Input", **style_label_header).grid(row=row1, column=0, **grid_label_header)
        row1+=1

        tk.Label(group, text="Sample table:", **style_label).grid(row=row1, column=0, **grid_label_col1st)
        self.sample_table_entry = tk.Entry(group, width=30, **style_entry)
        self.sample_table_entry.grid(row=row1, column=1, sticky='w', columnspan=2)
        tk.Button(group, text="Browse", command=lambda: select_file(self.sample_table_entry), **style_button_browse).grid(row=row1, column=3, sticky='w', padx=(10,0))

        self.set_option_default()

        tk.Button(group, text="RUN", command=self.run_var, **style_button_big).grid(row=100, column=0, pady=(20, 10), padx=(50,0), sticky='w', columnspan=2)
        tk.Button(group, text="RESET", command=self.reset_var, **style_button_big).grid(row=100, column=2, pady=(20, 10), sticky='w', columnspan=2)

        group.pack(anchor='center', expand=True)

    def set_option_default(self):
        self.detect_mode_combobox.set("germline")
        self.seq_method.set("wgs")
        self.file_type.set("fastq")
        self.vartype_variable[0].set(True)
        self.vartype_variable[1].set(False)
        self.vartype_variable[2].set(False)      

    def run_var(self):

        def run_var_shell():
            path = option["outdir"] + "/variant_calling.log"
            with open(path, "w") as log_file:
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while True:
                    output = proc.stdout.readline().decode().strip()
                    if not output and proc.poll() is not None:
                        break
                    output_box.config(state="normal")
                    output_box.insert(tk.END, output + "\n")
                    output_box.config(state="disabled")
                    output_box.see(tk.END)
                    log_file.write(output + "\n")
            log_file.close()

            proc.communicate()
            if (proc.returncode != 0):
                tk.messagebox.showinfo("INFO", "Execution failed. See log in output window or output directory.")
            else:
                tk.messagebox.showinfo("INFO", "Execution completed!")

        option = {
            "dect_mode": self.detect_mode_combobox.get(),
            "seq_method": self.seq_method.get(),
            "file_type": self.file_type.get(),
            
            "sample_table": self.sample_table_entry.get(),
            "region": self.entry_region.get(),
            "outdir": self.entry_outdir.get(),

            "short": self.vartype_variable[0].get(),
            "sv": self.vartype_variable[1].get(),
            "cnv": self.vartype_variable[2].get(),
        }
        var_type = ""
        if option["short"]:
            var_type+="snp,indel,"
        if option["sv"]:
            var_type+="sv,"
        if option["cnv"]:
            var_type+="cnv,"
        if not option["sample_table"]:
            tk.messagebox.showerror("ERROR", "No sample information table!")
        if option["seq_method"] == "wes" and not option["region"]:
            tk.messagebox.showerror("ERROR", "Region .bed file is required when analyzing wes or target data!")

        optional = {
            "--remove-dup": advance_option["rmdup"],
            "--soft-clip": advance_option["softclip"],
            "--force-rg": advance_option["forceRG"],
            "--annot": advance_option["default"],
            "--symbol": advance_option["symbol"],
            "--exist": advance_option["exist"],
            "--pathogenicity": advance_option["pathogenicity"],
            "--frequency": advance_option["frequency"],
            "--clinvar": advance_option["clinvar"],
            "--report": advance_option["report"]
        }
        dnapipe_script = dnapipe_dir + "/dnapipe"
        command = ["bash", dnapipe_script, "var", "--mode", option["dect_mode"], 
                "--seq-method", option["seq_method"], "--file-type", option["file_type"], 
                "--variant", var_type, "-o", option["outdir"], "--sample-info", option["sample_table"],
                "-t", advance_option["threads"], "-m", advance_option["memory"], "--mapq", advance_option["mapq"],
                "--short", advance_option["short_tools"], "--sv", advance_option["sv_tools"], "--cnv", advance_option["cnv_tools"]]
        command.extend([key for key, value in optional.items() if value])
        if option["region"]:
            command.extend(["--region", option['region']])
        if advance_option["bin"]:
            command.extend(["--bin", advance_option['bin']])

        t = threading.Thread(target=run_var_shell)
        t.start()

    def reset_var(self):
        entry_list = [self.entry_region, self.entry_outdir, self.sample_table_entry]
        for entry in entry_list:
            entry.delete(0,tk.END)
        self.set_option_default()
        
class visual_page(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        
        self.configure(bg="white")
        self.vcf_list = ["SNV:", "INDEL:", "SV:", "CNV:"]
        self.visual_list = ["Figure:", "","Analysis:"]
        self.plot_list = ["Category", "Length", "Pathogenicity*", "Circos"]
        self.analysis_list = ["GO*", "KEGG*", "PPI*"]
        self.plot_var = {}
        self.analysis_var = {}
        self.visual_group = tk.LabelFrame(self, text="",borderwidth=0, **style_group)
        self.create_widgets()

    def create_widgets(self):

        pos_label_header1 = {"sticky":'w', "columnspan": 2, "pady": 10}
        pos_label_col1st = {"sticky":'w', "pady": 2, "padx": (30,0)}

        tk.Label(self.visual_group, text="VCF file Input", **style_label_header).grid(row=0, column=0, **pos_label_header1)
        self.entries = []
        for index in range(len(self.vcf_list)):
            tk.Label(self.visual_group, text=self.vcf_list[index], **style_label).grid(row=index + 1, column=0, **pos_label_col1st)
            entry = tk.Entry(self.visual_group, width=30, **style_entry)
            entry.grid(row=index + 1, column=1, sticky='w', columnspan=2, padx=(0,10), pady=2)
            self.entries.append(entry)
            tk.Button(self.visual_group, text="Browse", command = lambda idx=index: select_file(self.entries[idx]), **style_button_browse).grid(row=index + 1, column=3, sticky='w', pady=2)
        
        style_checkbox2 = { "bg": "white", "highlightthickness":0 , "anchor": "w", "width": 12 }
        tk.Label(self.visual_group, text="Visualize about", **style_label_header).grid(row=5, column=0, **pos_label_header1)
        for index in range(len(self.visual_list)):
            tk.Label(self.visual_group, text=self.visual_list[index], **style_label).grid(row=index + 6, column=0, **pos_label_col1st)
        for i in range(len(self.plot_list)):
            self.plot_var[i] = tk.BooleanVar()
            if i != (len(self.plot_list)-1):
                tk.Checkbutton(self.visual_group, text=self.plot_list[i], variable=self.plot_var[i], **style_checkbox2).grid(row=6, column=i%3+1, padx=(0,10), sticky='w', pady=3)
        tk.Checkbutton(self.visual_group, text="Circos", variable=self.plot_var[len(self.plot_list)-1], command=self.circos_click, **style_checkbox2).grid(row=7, column=1, padx=(0,10), sticky='w', pady=3)
        self.bin_entry = tk.Entry(self.visual_group, width=10, **style_entry)
        self.bin_entry.insert(0, "5000000")
        
        for i in range(len(self.analysis_list)):
            self.analysis_var[i] = tk.BooleanVar()
            tk.Checkbutton(self.visual_group, text=self.analysis_list[i], variable=self.analysis_var[i], **style_checkbox2, command=self.analysis_click).grid(row=8, column=i+1, padx=(0,10), sticky='w', pady=3)
        self.entry_cadd = tk.Entry(self.visual_group, width=5, **style_entry)
        self.entry_cadd.insert(0, "10")

        tk.Label(self.visual_group, text="Output", **style_label_header).grid(row=15, column=0, **pos_label_header1)
        tk.Label(self.visual_group, text="Directory:", **style_label).grid(row=16, column=0, **pos_label_col1st)
        entry = tk.Entry(self.visual_group, width=30, **style_entry)
        entry.grid(row=16, column=1, sticky='w', columnspan=2, padx=(0,10))
        self.entries.append(entry)
        tk.Button(self.visual_group, text="Browse", command=lambda: select_folder(entry), **style_button_browse).grid(row=16, column=3, sticky='w')

        tk.Button(self.visual_group, text="Run", command=self.run_visual, **style_button_big).grid(row=99, column=0, sticky='nw', padx=(30, 0), pady=15, columnspan=2)
        tk.Button(self.visual_group, text="Reset", command=self.reset, **style_button_big).grid(row=99, column=2, sticky='nw', padx=(30, 0),  pady=15, columnspan=2)
 
        self.visual_group.pack(anchor='center', expand=True)
    
    def circos_click(self):
        if self.plot_var[len(self.plot_list)-1].get():
            self.label = tk.Label(self.visual_group, text="Bin size (bp):", **style_label)
            self.label.grid(row=7, column=2, padx=(0,10), sticky='w')
            self.bin_entry.grid(row=7, column=3, sticky='w')
        else:
            self.label.grid_forget()
            self.bin_entry.grid_forget()
    
    def analysis_click(self):
        if any(item for item in self.analysis_var):
            tk.Label(self.visual_group, text="CADD threshold:", **style_label).grid(row=9, column=1, sticky='w', padx=(0, 10), pady=2)
            self.entry_cadd.grid(row=9, column=2, sticky='w', padx=(0,10), pady=2)

    def run_visual(self):

        def run_visual_shell():
            path = self.input["out_dir"] + "/visual.log"
            with open(path, "w") as log_file:
                output_box.config(state="normal")
                output_box.insert(tk.END, "Plotting..." + "\n")
                output_box.config(state="disabled")
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                while True:
                    output = proc.stdout.readline().decode().strip()
                    if not output and proc.poll() is not None:
                        break
                    output_box.config(state="normal")
                    output_box.insert(tk.END, output + "\n")
                    output_box.config(state="disabled")
                    output_box.see(tk.END)
                    log_file.write(output + "\n")
            log_file.close()

            proc.communicate()
            tk.messagebox.showinfo("INFO", "Execution completed. Please check the results in the output directory and log information in the output window.")

        self.input = {
            "SNP": self.entries[0].get(),
            "INDEL": self.entries[1].get(),
            "SV": self.entries[2].get(),
            "CNV": self.entries[3].get(),
            
            "bin": self.bin_entry.get(),
            "cadd":self.entry_cadd.get(),
            "out_dir": self.entries[4].get()
        }

        plot_type = {
            "--category": self.plot_var[0].get(),
            "--length": self.plot_var[1].get(),
            "--pathogenicity": self.plot_var[2].get(),
            "--circos": self.plot_var[3].get(),
            "--go":self.analysis_var[0].get(),
            "--kegg":self.analysis_var[1].get(),
            "--ppi":self.analysis_var[2].get(),
        }
        
        if not( self.input["SNP"] or self.input["INDEL"] or self.input["CNV"] or self.input["SV"] ):
            tk.messagebox.showerror("ERROR", "No vcf input!") 
        elif not( self.input["out_dir"] ):
            tk.messagebox.showerror("ERROR", "Output directory not chosen!")
        elif self.plot_var[3].get() and not all(self.input.get(key) for key in ["SNP", "INDEL", "SV"]):
            tk.messagebox.showerror("ERROR", "Circos plot require at least SNP, INDEL and SV vcf files!")
        elif not self.input["bin"] and self.input["cadd"]:
            tk.messagebox.showerror("ERROR", "Empty value in bin size or cadd socre threshold!")
        else:               
            for key in ["SNP", "INDEL", "SV", "CNV"]:
                if not self.input[key]:
                    self.input[key] = "NA"
            dnapipe_script = dnapipe_dir + "/dnapipe"
            vcfs = ["--snp", self.input["SNP"], "--indel", self.input["INDEL"], "--sv", self.input["SV"], "--cnv", self.input["CNV"]]
            options = ["--bin", self.input["bin"], "--cadd_score", self.input["cadd"], "--outdir", self.input["out_dir"]]
            options.extend([key for key, value in plot_type.items() if value])
            command = ["bash", dnapipe_script, "plot"]
            command.extend(vcfs)
            command.extend(options)
            t = threading.Thread(target=run_visual_shell)
            t.start()
    
    def reset(self):
        for i in range(len(self.plot_list)):
            self.plot_var[i].set("False")
        for i in range(len(self.analysis_list)):
            self.analysis_var[i].set("False")
        for entry in self.entries:
            entry.delete(0,tk.END)
        self.bin_entry.delete(0,tk.END)
        self.bin_entry.insert(0, "5000000")
        self.entries[4].insert(0, "10")



class App(tk.Tk):

    def __init__(self):
        super().__init__()
        self.set_window()
        self.set_main()
        self.set_output()
        
    def set_window(self):
        s = ttk.Style()
        s.theme_use('default')
        s.configure('TNotebook', background="white")
        s.configure('TNotebook.Tab', background='whitesmoke')
        s.map("TNotebook.Tab", background= [("selected", "white")])
        s.configure('TCombobox', selectbackground="transparent", selectforeground="black", background="skyblue", arrowcolor="white")
        s.map('TCombobox', fieldbackground=[('readonly', 'whitesmoke')])
        self.title("hDNApipe")
        self.geometry("1000x500")

    def set_output(self):
        group1 = tk.LabelFrame(self, text="Output Window", bg="white", padx=5, pady=10, labelanchor="nw")
        frame_output = tk.Frame(group1, bg="white")
        frame_output.pack(side="right", fill="both", expand=True)
        scrollbar = tk.Scrollbar(frame_output, orient='vertical', bg="skyblue", troughcolor="white")
        scrollbar.pack(side=tk.RIGHT, fill='y')
        global output_box
        output_box = tk.Text(frame_output, wrap=tk.WORD, yscrollcommand=scrollbar.set, state="disabled")
        output_box.pack(fill="both", expand=True)
        scrollbar.config(command=output_box.yview)
        group1.pack(side="right", fill="both", expand=True)

    def set_main(self):
        frame_main = tk.Frame(self, bg="white", height=500)
        frame_main.pack(side="left", fill="both", expand=True)
        self.notebook = ttk.Notebook(frame_main, width=520)
        page_help = init_page(self.notebook)
        page_var = var_page(self.notebook)
        page_adv = advance_page(self.notebook)
        page_visual = visual_page(self.notebook)
        self.notebook.add(page_help, text="Initialization")
        self.notebook.add(page_var, text="Basic options")
        self.notebook.add(page_adv, text="Advanced options")
        self.notebook.add(page_visual, text="VCF visulization")
        self.notebook.pack(expand=True, fill=tk.BOTH)


if __name__=="__main__":

    app = App()
    app.mainloop()
