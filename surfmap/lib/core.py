#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Tuple, Union
import re
import json

from surfmap import PATH_MSMS, __COPYRIGHT_FULL__
from surfmap.lib.logs import get_logger
from surfmap.lib.parameters import Parameters
from surfmap.lib.utils import JunkFilePath
from surfmap.tools.SurfmapTools import run_particles_mapping
from surfmap.tools.compute_shell import run as run_compute_shell
from surfmap.tools.compute_electrostatics import run as run_compute_electrostatics
from surfmap.tools import Structure


logger = get_logger(name=__name__)


def generate_mab_tagging_resfile(pdb_file: Union[str, Path], mab_tagging: str, outdir: Union[str, Path]) -> str:
    """
    Parses a PDB file and generates a SURFMAP-compatible residue highlight file
    for the specified IMGT CDR regions. Supports specific chains (e.g., CDR3H, CDR1L)
    or all chains.
    """
    logger.info(f"MAb Tagging: Initializing extraction for tags '{mab_tagging}' using IMGT numbering...")
    dPDB = Structure.parsePDBMultiChains(str(pdb_file))
    
    imgt_regions = {
        "CDR1": range(27, 39),
        "CDR2": range(56, 66),
        "CDR3": range(105, 118)
    }
    
    tags = [t.strip().upper() for t in mab_tagging.split(':')]
    target_map = {} 
    
    for tag in tags:
        if tag.startswith("CDR"):
            region = tag[:4] 
            chain_id = tag[4:] 
            
            if region in imgt_regions:
                if chain_id:
                    for r in imgt_regions[region]:
                        target_map[(chain_id, r)] = tag
                    logger.debug(f"MAb Tagging: Mapped tag '{tag}' to region {region} in chain '{chain_id}'")
                else:
                    for chain in dPDB["chains"]:
                        for r in imgt_regions[region]:
                            target_map[(chain, r)] = tag
                    logger.debug(f"MAb Tagging: Mapped tag '{tag}' to region {region} across ALL chains")
            else:
                logger.warning(f"MAb Tagging: Unrecognized region '{region}' in tag '{tag}'. Ignored.")
        else:
            logger.warning(f"MAb Tagging: Unrecognized tag format '{tag}'. Expected format like CDR3H. Ignored.")
             
    resfile_path = Path(outdir) / f"{Path(pdb_file).stem}_mab_tagging.txt"
    found_count = 0
    
    with open(resfile_path, "w") as f:
        for chain in dPDB["chains"]:
            logger.info(f"MAb Tagging: Scanning chain '{chain}'...")
            chain_match_count = 0
            
            for res in dPDB[chain]["reslist"]:
                match = re.search(r'\d+', res)
                if match:
                    resid_int = int(match.group())
                    if (chain, resid_int) in target_map:
                        tag = target_map[(chain, resid_int)]
                        resname = dPDB[chain][res]["resname"]
                        f.write(f"{chain}\t{res}\t{resname}\t{tag}\n")
                        chain_match_count += 1
                        found_count += 1
                        logger.trace(f"MAb Tagging: Found match -> Chain {chain}, Resid {res}, Resname {resname}, Tag {tag}")
                        
            logger.info(f"MAb Tagging: Found {chain_match_count} matching residues in chain '{chain}'.")

    logger.info(f"MAb Tagging: Extraction complete. Found a total of {found_count} matching residues.")
    if found_count == 0:
        logger.warning("MAb Tagging: No residues matched the requested tags! Please verify that the PDB uses the IMGT numbering scheme and contains the requested chains.")
        
    return str(resfile_path)


def compute_coords_list(params: Parameters, coords_file: str, property: str) -> Tuple[int, str]:
    cmd = ["Rscript", params.coords_script, "-f", str(coords_file), "-s", str(params.cellsize), "-P", params.proj, "-o", params.outdir]
    logger.debug(f"Running the command: {' '.join(cmd)}")
    proc_status = subprocess.call(cmd)
    
    if proc_status != 0:
        logger.error(f"Error occured during the computation of spherical coordinates of particles, the process will stop.")
        exit(1)

    out_file = str(Path(params.outdir) / "coord_lists" / f"{Path(params.pdbarg).stem}_{property}_coord_list.txt")
    return proc_status, out_file


def generate_json_matrix(txt_path: str, pdb_path: str):
    """
    Parses a generated text matrix and maps contributing residues (if any) to their sequential numbering
    1..N inside their respective chains, returning a JSON equivalent of the grid.
    """
    if not Path(txt_path).exists():
        return
    
    json_path = txt_path.replace(".txt", ".json")
    out_data = []
    
    # Safely convert R's NA/Inf strings to valid Python floats or None
    def _safe_float(val_str):
        v = val_str.strip().upper()
        if v in ('INF', '-INF', 'NA', 'NAN', ''):
            return None
        try:
            return float(val_str.strip())
        except ValueError:
            return None
    
    # Build mapping from chain and PDB resid (including insertion codes) -> sequential number
    try:
        dPDB = Structure.parsePDBMultiChains(str(pdb_path))
        seq_map = {}
        for chain in dPDB.get("chains", []):
            seq_counter = 1
            for res in dPDB[chain].get("reslist", []):
                seq_map[(chain, res)] = seq_counter
                seq_counter += 1
    except Exception as e:
        logger.warning(f"Could not parse PDB for JSON mapping: {e}")
        seq_map = {}

    try:
        with open(txt_path, "r") as f:
            lines = f.readlines()
            if not lines:
                return
                
            for line in lines[1:]: # skip header
                parts = line.strip('\n').split('\t')
                if len(parts) < 4:
                    continue
                absc, ord_val, val, residues = parts[0], parts[1], parts[2], parts[3]
                
                res_list = []
                if residues != "NA" and residues.strip() != "":
                    for r in residues.split(","):
                        r = r.strip()
                        if not r: continue
                        r_parts = r.split("_")
                        if len(r_parts) >= 3:
                            resname = r_parts[0]
                            resnum = r_parts[1]
                            chain_id = r_parts[2]
                            seq_num = seq_map.get((chain_id, resnum), None)
                            
                            res_list.append({
                                "resname": resname,
                                "pdb_resnum": resnum,
                                "chain": chain_id,
                                "seq_num": seq_num
                            })
                
                out_data.append({
                    "absc": _safe_float(absc),
                    "ord": _safe_float(ord_val),
                    "value": _safe_float(val),
                    "residues": res_list
                })
                
        with open(json_path, "w") as jf:
            json.dump(out_data, jf, indent=2)
    except Exception as e:
        logger.warning(f"Failed to generate JSON for {txt_path}: {e}")


def compute_matrix(params: Parameters, coords_file: str, property: str, suffix="_coord_list.txt") -> Tuple[int, str, str]:
    cmd = ["Rscript", params.matrix_script, "-i", coords_file, "-s", str(params.cellsize), "-P", str(params.proj), "-o", str(params.outdir), "--suffix", suffix, "--discrete"]

    if property != "binding_sites":
        if not params.nosmooth:
            del cmd[-1]
        else:
            cmd[-1] = "--nosmooth"

    logger.debug(f"Running the command: {' '.join(cmd)}")
    proc_status = subprocess.call(cmd)

    if proc_status != 0:
        logger.error(f"Error occured during smoothing of the raw matrix, the process will stop.")
        exit(1)

    named_property = "bfactor" if property == "binding_sites" else property
    out_matrix_smoothed = str(Path(params.outdir) / "smoothed_matrices" / f"{Path(params.pdbarg).stem}_{named_property}_smoothed_matrix.txt")
    out_matrix = str(Path(params.outdir) / "matrices" / f"{Path(params.pdbarg).stem}_{named_property}_matrix.txt")

    # Generate JSON structured equivalents of both matrices
    if params.pdbarg and Path(params.pdbarg).exists():
        generate_json_matrix(out_matrix, params.pdbarg)
        generate_json_matrix(out_matrix_smoothed, params.pdbarg)

    return proc_status, out_matrix, out_matrix_smoothed        


def compute_map(params: Parameters, matrix_file: str, property: str, reslist: str=None, suffix="_smoothed_matrix.txt") -> int:
    out_pdf = str(Path(params.outdir) / "maps" / f"{params.pdb_id}_{property}_map.pdf")
    out_png = None

    if property == "binding_sites":
        scale_opt = "--discrete"
    elif property == "circular_variance_atom":
        scale_opt = "--circular_variance"
    else:
        scale_opt = "--" + property

    cmd = ["Rscript", params.map_script, "-i", str(matrix_file), scale_opt,  "-s", str(params.cellsize), "-p", params.pdb_id, "-P", params.proj, "-o", str(params.outdir), "--suffix", suffix]
    
    if params.coordstomap:
        cmd += ["-c", params.coordstomap]
    if reslist:
        cmd += ["-l", reslist]
    if params.png:
        cmd.append("--png")
        out_png = out_pdf.replace(".pdf", ".png")

    if getattr(params, 'margin_scale', 1.0) != 1.0:
        cmd += ["--margin_scale", str(getattr(params, 'margin_scale'))]
    if getattr(params, "no_scale_bar", False):
        cmd.append("--no_scale_bar")
    if params.elec_max_value is not None:
        cmd += ["--elec_max_value", str(params.elec_max_value)]
    if params.bfactor_min_value is not None:
        cmd += ["--bfactor_min_value", str(params.bfactor_min_value)]
    if params.bfactor_max_value is not None:
        cmd += ["--bfactor_max_value", str(params.bfactor_max_value)]

    logger.debug(f"Running the command: {' '.join(cmd)}")
    proc_status = subprocess.call(cmd)

    if proc_status != 0:
        logger.error(f"Error occured during computing the map, the process will stop.")
        exit(1)

    return proc_status, out_png, out_pdf


def surfmap_from_pdb(params: Parameters, with_copyright: bool=True):
    if with_copyright:
        print(__COPYRIGHT_FULL__)
    
    junk_optional = JunkFilePath(elements=[Path(params.outdir) / "shells", Path(params.outdir) / "tmp-elec"])
    shell = None

    if hasattr(params, 'mab_tagging') and params.mab_tagging:
        mab_resfile = generate_mab_tagging_resfile(params.pdbarg, params.mab_tagging, params.outdir)
        if not params.resfile:
            params.resfile = mab_resfile
        else:
            with open(mab_resfile, 'r') as f_in:
                mab_data = f_in.read()
            with open(params.resfile, 'a') as f_out:
                f_out.write("\n" + mab_data)
        junk_optional.add(element=[mab_resfile])

    listtomap = ["kyte_doolittle", "stickiness", "wimley_white", "circular_variance"] if params.ppttomap == "all" else [params.ppttomap]     
    
    for tomap in listtomap:
        logger.info(msg=f"Surface mapping of the {tomap} property".upper())

        step_index = 1
        if shell is None:
            outdir_shell = Path(params.outdir) / "shells"
            extra_radius = params.rad
            logger.info(msg=f"Step {step_index}: computing a shell around the protein surface")
            csv_coords, shell = run_compute_shell(pdb_filename=params.pdbarg, out_dir=outdir_shell, extra_radius=extra_radius)
        else:
            logger.info(msg=f"Step {step_index}: shell already exists. Skipping this step")

        if params.ppttomap == "electrostatics":
            step_index += 1
            outdir_elec = Path(params.outdir) / "tmp-elec"
            logger.info(msg=f"Step {step_index}: computing electrostatics potential")
            shell = run_compute_electrostatics(pdb_filename=params.pdbarg, csv_filename=csv_coords, force_field=params.force_field, pqr_filename=params.pqr, out_dir=outdir_elec)

        step_index += 1
        logger.info(msg=f"Step {step_index}: computing the property values and/or assign it to the shell particles")
        property = "bfactor" if tomap == "binding_sites" else tomap
        reslist, partlist_outfile = run_particles_mapping(shell=shell, pdb=params.pdbarg, tomap=property, outdir=params.outdir, res=params.resfile)
        junk_optional.add(element=[reslist, partlist_outfile])

        step_index += 1
        logger.info(msg=f"Step {step_index}: computing the 2D {params.proj} projection coordinates of each shell particle")
        _, coordfile = compute_coords_list(params=params, coords_file=partlist_outfile, property=property)
        junk_optional.add(element=[coordfile, Path(coordfile).parent])
        
        step_index += 1
        logger.info(msg=f"Step {step_index}: dividing the 2D projection into {int(360/params.cellsize)}x{int(180/params.cellsize)} cells and smoothing the values")
        _, matrix_file, smoothed_matrix_file = compute_matrix(params=params, coords_file=coordfile, property=tomap)
        junk_optional.add(element=[matrix_file, Path(matrix_file).parent])

        step_index += 1
        logger.info(msg=f"Step {step_index}: computing the 2D map")
        _, png_filename, pdf_filename = compute_map(params=params, matrix_file=smoothed_matrix_file, property=tomap, reslist=reslist)

    if not params.keep:
        junk_optional.empty()
        
    params.write_parameters(filename="parameters.log")
    print()


def surfmap_from_matrix(params: Parameters, with_copyright: bool=True):
    if with_copyright:
        print(__COPYRIGHT_FULL__)

    if params.ppttomap == 'all':
        logger.error("Error: the property to map cannot be set to 'all' when computing a map from a matrix file.\n")
        exit()

    matrices_outdir = Path(params.outdir) / "matrices"
    matrices_outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(params.mat, matrices_outdir / Path(params.mat).name)

    matf = Path(params.outdir) / "matrices" / Path(params.mat).name
    
    logger.info(msg=f"Computing the 2D map")
    _, png_filename, pdf_filename = compute_map(params=params, matrix_file=matf, property=params.ppttomap)
    
    params.write_parameters(filename="parameters.log")
