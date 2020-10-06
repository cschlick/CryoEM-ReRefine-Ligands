
# phenix.ligands_cc

> A tool for CCTBX/Phenix to extract ligands and calculate map-model cross correlation.
> It will optionally write out density maps and models for each ligand

Usage: phenix.ligands_cc <model_file> <map_file> resolution= <resolution>  

Optional command line arguments:

 write_ligand_files = False  
      .type = bool  
      .help = Whether to write the map and model for each ligand  

 output_directory = None  
      .type = str  
      .help = Optional directory to write ligand maps and models  

 include_only_resnames = None  
      .type = str  
      .multiple = True  
      .help = Only extract ligands with this residue name  

 include_resnames = None  
      .type = str  
      .multiple = True  
      .help = Do extract ligands with this residue name  

 exclude_resnames = None  
      .type = str  
      .multiple = True  
      .help = Do NOT extract ligands with this residue name  

 include_residue_classes = None  
      .type = str  
      .multiple = True  
      .help = Do extract ligands of this residue class  
 
 exclude_residue_classes = None  
      .type = str  
      .multiple = True  
      .help = Do NOT extract ligands of this residue class  




