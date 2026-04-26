import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, AllEnzymes

def load_fasta(file_path):
    """Wczytuje pierwszą sekwencję z pliku FASTA."""
    for record in SeqIO.parse(file_path, "fasta"):
        return record.seq
    return None

def perform_restriction_analysis(seq):
    """
    Wykonuje analizę restrykcyjną z użyciem dostępnych enzymów.
    Zwraca unikalne długości fragmentów DNA (bez powtórzeń).
    """
    enzymes = RestrictionBatch(AllEnzymes)
    analysis = Analysis(enzymes, seq)
    results = analysis.full()

    fragment_results = {}

    for enzyme, cuts in results.items():
        if cuts:
            cuts_sorted = sorted(cuts)
            fragments = []

            prev_cut = 0
            for cut in cuts_sorted:
                fragments.append(cut - prev_cut)
                prev_cut = cut

            fragments.append(len(seq) - prev_cut)

            # Usunięcie duplikatów i sortowanie
            unique_fragments = sorted(set(fragments))

            fragment_results[str(enzyme)] = unique_fragments

    return fragment_results

def save_results(results, output_path):
    """Zapisuje wyniki do pliku tekstowego."""
    with open(output_path, "w") as f:
        for enzyme, fragments in results.items():
            f.write(f"{enzyme}:\n")
            f.write("Fragment sizes (bp): " + ", ".join(map(str, fragments)) + "\n\n")

def run_analysis():
    fasta_path = filedialog.askopenfilename(
        title="Wybierz plik FASTA",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )

    if not fasta_path:
        return

    seq = load_fasta(fasta_path)
    if seq is None:
        messagebox.showerror("Błąd", "Nie udało się wczytać sekwencji.")
        return

    results = perform_restriction_analysis(seq)

    save_path = filedialog.asksaveasfilename(
        title="Zapisz wyniki",
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt")]
    )

    if not save_path:
        return

    save_results(results, save_path)

    messagebox.showinfo("Sukces", "Analiza zakończona i zapisana do pliku.")

# GUI
root = tk.Tk()
root.title("Analiza restrykcyjna DNA")

frame = tk.Frame(root, padx=20, pady=20)
frame.pack()

btn = tk.Button(frame, text="Wybierz plik FASTA i analizuj", command=run_analysis)
btn.pack()

root.mainloop()
