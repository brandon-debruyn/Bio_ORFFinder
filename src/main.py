# Author: Brandon Lee De Bruyn
# Date: 12/04/2024
# Description: A basic TKinter GUI implementation for the ORF Finder application. 
# Note: Allows the user to enter a DNA sequence manually or upload a file, and then find the open reading frames in the sequence. 
# There is CLI output as well as "in-window" output in the GUI. (both incase the user wants to use the GUI or CLI for the application)

import os
import tkinter as tk
from tkinter import scrolledtext, filedialog, ttk

from orf_finder import ORFFinder
from helpers import StartCodon

class ORFApp:
    def __init__(self, master):
        # set the master window
        self.master = master
        master.title("ORF Finder") # set the title of the window

        # create a label widget to display upload text
        self.label = tk.Label(master, text="Enter Nucleotide Sequence Manually or Upload a File:")
        self.label.pack()

        # create a text widget to allow the user to enter a nucleotide sequence
        self.nucleotide_sequence = tk.Text(master, height=15, width=150) 
        self.nucleotide_sequence.pack()

        # create a button widget to allow the user to upload a file
        self.upload_button = tk.Button(master, text="Upload Nucleotide Sequence File", command=self.upload_file) 
        self.upload_button.pack()

        # label for minimum ORF length 
        self.min_orf_label = tk.Label(master, text="Select Minimum ORF Length (default 30):")
        self.min_orf_label.pack()

        # dropdown for selecting minimum ORF length
        self.min_orf_length = tk.StringVar()
        self.min_orf_combo = ttk.Combobox(master, textvariable=self.min_orf_length)
        self.min_orf_combo['values'] = (30, 75, 150)
        self.min_orf_combo['state'] = 'readonly'  # readonly combo
        self.min_orf_combo.pack()
        self.min_orf_combo.current(0)  # default to 30

        # label for start codon selection
        self.start_codon_label = tk.Label(master, text="Start Codon:")
        self.start_codon_label.pack()

        # variable to track the start codon selection
        self.start_codon_var = tk.StringVar(value=StartCodon.AUTOMATIC.value)

        # iteratively create radiobuttons for selecting the start codon
        for codon in StartCodon:
            rb = tk.Radiobutton(master, text=codon.value, variable=self.start_codon_var, value=codon.value)
            rb.pack()

        # create a button widget to allow the user to find ORFs in the nucleotide sequence
        self.find_orf_button = tk.Button(master, text="Find ORFs", command=self.find_orfs) 
        self.find_orf_button.pack()

        # create a scrolled text widget to display the ORFs found in the DNA sequence
        self.orf_display = scrolledtext.ScrolledText(master, height=15, width=150) 
        self.orf_display.pack()
    
    def upload_file(self):

        # set the allowed file extensions
        allowed_extensions = ['.txt', '.fasta', '.fa']
    
        # open a file dialog to allow the user to select a file
        file_path = filedialog.askopenfilename()
        if file_path:
            # extract and check if the file type is allowed
            _, file_extension = os.path.splitext(file_path)
            if file_extension.lower() not in allowed_extensions:

                # print an error message to the console and exit the function if the file type is not allowed
                print(f"Unsupported file type: {file_extension}. Please upload a .txt, .fasta, or .fa file.") 
                return  
            
            # read the file and insert the nucleotide sequence into the text widget
            try:
                with open(file_path, 'r') as file:
                    nucleotide_sequence = file.read().strip()
                    self.nucleotide_sequence.delete("1.0", tk.END)
                    self.nucleotide_sequence.insert(tk.END, nucleotide_sequence)
                    self.orf_display.delete("1.0", tk.END)  # Clear the ORF display area
            except Exception as e:  # catch any exceptions that occur when reading the file
                print(f"Error reading the file: {e}")
                self.nucleotide_sequence.delete("1.0", tk.END)  
                self.nucleotide_sequence.insert(tk.END, e.__str__()) # insert the error message part of the error object into the text widget  

    def find_orfs(self):

        # clear the ORF display
        self.orf_display.delete("1.0", tk.END)
        
        # get the nucleotide sequence, minimum ORF length, and start codon selection
        nucleotide_sequence = self.nucleotide_sequence.get("1.0", tk.END).strip()

        min_orf_length = int(self.min_orf_length.get())
        selected_start_codon = StartCodon(self.start_codon_var.get())

        # create an ORFFinder object and find the ORFs
        orf_finder = ORFFinder(nucleotide_sequence, min_orf_length, selected_start_codon)
        orf_data, display_data = orf_finder.find_orfs()
            
        # insert the ORF data
        self.orf_display.insert(tk.END, display_data)
        
# run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = ORFApp(root)
    root.mainloop()