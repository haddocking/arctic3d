# RODIGY = "/trinity/login/jlopez/software/prodigy-cryst/interface_classifier.py" # Location of Prodigy Crystal


# def run_prodigy(pdb_f, chain_i, chain_j):
#     """Wrapper for running PRODIGY-CRYSTAL."""

#     cmd = f"python {PRODIGY} {pdb_f} --selection {chain_i} {chain_j}"

#     p = subprocess.Popen(
#         cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
#     )

#     out, _ = p.communicate()
#     try:
#         interface_class = out.decode("utf-8").split(os.linesep)[-2].split()[2]
#     except Exception as e:
#         logging.debug(f"{pdb_f}, {e}")
#         interface_class = "NONE"

#     return interface_class


# def validate_interaction(prot_a, prot_b, pre_sel_pdbs, prodigy_results_pdb, data_folder):
#    """Check if the interaction between two proteins is biological."""
#    biological_pdbs = []
#    # Interface data Protein A
#    url_a = f"{INTERFACE_URL}/{prot_a}"
#    try:
#        interface_data_a = make_request(url_a, None)
#    except Exception as e:
#        logging.debug(f"Could not make InterfaceResidues request for {prot_a}, {e}")
#        return None
#
#    # Interface data Protein B
#    url_b = f"{INTERFACE_URL}/{prot_b}"
#    try:
#        interface_data_b = make_request(url_b, None)
#    except Exception as e:
#        logging.debug(f"Could not make InterfaceResidues request for {prot_b}, {e}")
#        return None
#
#    # Check if has been processed totally or partially before
#    protein_pair = f"{prot_a}-{prot_b}"
#    if protein_pair not in prodigy_results_pdb.keys(): # The pair has not been processed before
#        prodigy_results_pdb[protein_pair] = {}
#    else:
#        # The pair has been processed before
#        not_redo_pdbs = []
#        for pre_sel_pdb in pre_sel_pdbs:
#            if pre_sel_pdb in prodigy_results_pdb[protein_pair].keys(): # The pdb for the pair has been processed before
#                not_redo_pdbs += [pre_sel_pdb]
#        # Avoid redoing already done pdbs.
#        for not_redo_pdb in not_redo_pdbs:
#            if not_redo_pdb in pre_sel_pdbs:
#                pre_sel_pdbs.remove(not_redo_pdb)
#
#    # Obtain a dictionary with pdb_id : chain_id for both proteins. *Jesus: I used a previous script for this
#    pdb_chain_a = dict(
#        [
#            (k["pdbId"], k["chainIds"])
#            for e in interface_data_a[prot_a]["data"]
#            if e["accession"] == prot_b
#            for j in e["residues"]
#            for k in j["interactingPDBEntries"]
#            if k["pdbId"] in pre_sel_pdbs
#        ]
#    )
#
#    pdb_chain_b = dict(
#        [
#            (k["pdbId"], k["chainIds"])
#            for e in interface_data_b[prot_b]["data"]
#            if e["accession"] == prot_a
#            for j in e["residues"]
#            for k in j["interactingPDBEntries"]
#            if k["pdbId"] in pre_sel_pdbs
#        ]
#    )
#
#    # Check that the pdbs are the same for both proteins
#    checked_pdb_chain_a = {}
#    checked_pdb_chain_b = {}
#
#    for pdb in pdb_chain_a.keys():
#        if pdb in pdb_chain_b.keys():
#            checked_pdb_chain_a[pdb] = pdb_chain_a[pdb]
#            checked_pdb_chain_b[pdb] = pdb_chain_b[pdb]
#    assert checked_pdb_chain_a.keys() == checked_pdb_chain_b.keys()
#
#    # Run Prodigy
#    for pdb in checked_pdb_chain_a:
#        renum_pdb_f = Path(f"{data_folder}/{pdb}_renum.pdb")
#        if renum_pdb_f.exists(): # Check if the pdb has already been downloaded and renumbered
#            pdb_f = renum_pdb_f
#        else: # Check if the pdb has already been downloaded
#            pdb_f = Path(f"{data_folder}/{pdb}.pdb")
#
#        if not pdb_f.exists(): # Download the pdb. *Jesus: this can be done by PDBrenum too
#            cmd = f'pdb_fetch {pdb} | grep "ATOM" | pdb_tidy'
#            p = subprocess.Popen(
#                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
#            )
#
#            out, err = p.communicate()
#
#            with open(pdb_f, "w") as fh:
#                fh.write(out.decode("utf-8"))
#
#        chain_a_list = checked_pdb_chain_a[pdb].split(",")
#        chain_b_list = checked_pdb_chain_b[pdb].split(",")
#        job_list = []
#        pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
#        for chain_a in chain_a_list:
#            for chain_b in chain_b_list:
#                if chain_a != chain_b:
#                    # Prodigy wont work if chains are the same..!
#                    job_list.append(
#                        (
#                            pdb_f,
#                            chain_a,
#                            chain_b,
#                            pool.apply_async(
#                                run_prodigy, args=(pdb_f, chain_a, chain_b)
#                            ),
#                        )
#                    )
#
#        for element in job_list:
#            pdb_f, chain_a, chain_b, proc = element
#            interface_class = proc.get()
#            PDB_DB.append((prot_a, prot_b, pdb_f, chain_a, chain_b, interface_class))
#            logging.info(
#                f"===> {prot_a}_{prot_b} = {pdb} {chain_a}-{chain_b} "
#                f"= {interface_class}"
#            )
#            if interface_class == "BIO":
#                # If any chain-chain interaction for this pdb is biological, the evaluation is biological
#                if pdb not in prodigy_results_pdb[protein_pair].keys():
#                   prodigy_results_pdb[protein_pair][pdb] = "biological"
#
#        pool.close()
#
#    for pdb in pre_sel_pdbs:
#        if pdb not in prodigy_results_pdb[protein_pair].keys(): # Those pdbs haven't been evaluated as biological before, so they are not
#            prodigy_results_pdb[protein_pair][pdb] = "NOT biological"
#
#    # Confirm biological pdbs
#    for processsed_pdb in prodigy_results_pdb[protein_pair].keys():
#        if prodigy_results_pdb[protein_pair][processsed_pdb] == "biological" and processsed_pdb in pre_sel_pdbs:
#            biological_pdbs += [processsed_pdb]
#
#    if len(biological_pdbs) == 0:
#        biological_pdbs = None
#
#    return biological_pdbs, prodigy_results_pdb
