#!/usr/bin/env python3
import os
import random
import argparse
import sys

# Constantes por defecto
DEFAULT_NUM_PROTEINS       = 100
DEFAULT_NUM_LIGANDS        = 100
DEFAULT_NUM_ATOMS_PROTEIN  = 1000
DEFAULT_NUM_ATOMS_LIGAND   = 100
DEFAULT_PROTEIN_DIR        = "data/proteins"
DEFAULT_LIGAND_DIR         = "data/ligands"

def print_help():
    help_message = f"""
Uso: generate_dataset.py [opciones]

Opciones:
  --num_proteins N         Número de proteínas a generar (default: {DEFAULT_NUM_PROTEINS})
  --num_ligands N          Número de ligandos a generar (default: {DEFAULT_NUM_LIGANDS})
  --num_atoms_protein N    Número de átomos por proteína (default: {DEFAULT_NUM_ATOMS_PROTEIN})
  --num_atoms_ligand N     Número de átomos por ligando (default: {DEFAULT_NUM_ATOMS_LIGAND})
  --protein_dir DIR        Directorio de salida para proteínas (default: "{DEFAULT_PROTEIN_DIR}")
  --ligand_dir DIR         Directorio de salida para ligandos (default: "{DEFAULT_LIGAND_DIR}")
  -h, --help               Mostrar este mensaje de ayuda y salir
"""
    print(help_message)
    sys.exit(0)

def generate_protein(filename, num_atoms):
    with open(filename, "w") as f:
        # Encabezado PDB simplificado
        f.write("HEADER    GENERATED PROTEIN\n")
        f.write("TITLE     Random Protein\n")
        for i in range(1, num_atoms + 1):
            # Genera coordenadas aleatorias entre -100.0 y 100.0
            x = random.uniform(-100.0, 100.0)
            y = random.uniform(-100.0, 100.0)
            z = random.uniform(-100.0, 100.0)
            # Construir la línea: las coordenadas se ubican en columnas fijas
            prefix = "ATOM" + " " * (30 - len("ATOM"))
            coord_str = f"{x:8.3f}{y:8.3f}{z:8.3f}"
            line = prefix + coord_str
            line = line.ljust(80)  # Completar a 80 caracteres
            f.write(line + "\n")
        f.write("TER\nEND\n")

def generate_ligand(filename, num_atoms):
    with open(filename, "w") as f:
        # Encabezado SDF simplificado
        f.write("Generated Ligand\n")
        f.write("Programmatically generated\n")
        f.write("Comment line\n")
        count_line = f"{num_atoms:3d}  0  0  0  0  0            999 V2000"
        f.write(count_line + "\n")
        for i in range(num_atoms):
            # Genera coordenadas aleatorias entre -50.0 y 50.0
            x = random.uniform(-50.0, 50.0)
            y = random.uniform(-50.0, 50.0)
            z = random.uniform(-50.0, 50.0)
            element = random.choice(["C", "N", "O", "H", "S"])
            # Formato SDF: cada campo de coordenada ocupa 10 caracteres
            line = f"{x:10.4f}{y:10.4f}{z:10.4f} {element:>3}"
            f.write(line + "\n")
        f.write("M  END\n")
        f.write("$$$$\n")

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--num_proteins', type=int, default=DEFAULT_NUM_PROTEINS,
                        help=f'Número de proteínas a generar (default: {DEFAULT_NUM_PROTEINS})')
    parser.add_argument('--num_ligands', type=int, default=DEFAULT_NUM_LIGANDS,
                        help=f'Número de ligandos a generar (default: {DEFAULT_NUM_LIGANDS})')
    parser.add_argument('--num_atoms_protein', type=int, default=DEFAULT_NUM_ATOMS_PROTEIN,
                        help=f'Número de átomos por proteína (default: {DEFAULT_NUM_ATOMS_PROTEIN})')
    parser.add_argument('--num_atoms_ligand', type=int, default=DEFAULT_NUM_ATOMS_LIGAND,
                        help=f'Número de átomos por ligando (default: {DEFAULT_NUM_ATOMS_LIGAND})')
    parser.add_argument('--protein_dir', type=str, default=DEFAULT_PROTEIN_DIR,
                        help=f'Directorio de salida para proteínas (default: "{DEFAULT_PROTEIN_DIR}")')
    parser.add_argument('--ligand_dir', type=str, default=DEFAULT_LIGAND_DIR,
                        help=f'Directorio de salida para ligandos (default: "{DEFAULT_LIGAND_DIR}")')
    parser.add_argument('-h', '--help', action='store_true',
                        help='Mostrar este mensaje de ayuda')
    
    args, unknown = parser.parse_known_args()
    
    if args.help:
        print_help()
    
    protein_dir = args.protein_dir
    ligand_dir = args.ligand_dir
    num_proteins = args.num_proteins
    num_ligands = args.num_ligands
    num_atoms_protein = args.num_atoms_protein
    num_atoms_ligand = args.num_atoms_ligand

    os.makedirs(protein_dir, exist_ok=True)
    os.makedirs(ligand_dir, exist_ok=True)

    print("Generando archivos de proteínas...")
    for i in range(1, num_proteins + 1):
        filename = os.path.join(protein_dir, f"protein_{i:03d}.pdb")
        generate_protein(filename, num_atoms_protein)
        
    print("Generando archivos de ligandos...")
    for i in range(1, num_ligands + 1):
        filename = os.path.join(ligand_dir, f"ligand_{i:03d}.sdf")
        generate_ligand(filename, num_atoms_ligand)
        
    print("Generación del dataset completada.")

if __name__ == "__main__":
    main()
