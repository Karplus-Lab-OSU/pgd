-- DATABASE conversion script.  This script converts data in the existing tables to the new table structure
--
-- NOTE THAT THIS SCRIPT OVERWRITES *ALL* DATA IN THE DJANGO MODELS FOR PGD_CORE!!
-- use caution when running this.

-- create mapping table for aa
DROP TABLE IF EXISTS aa_codes;
CREATE TABLE `aa_codes` (
    `key` varchar(4) NOT NULL PRIMARY KEY,
    `value` varchar(4) NOT NULL
);

-- enter aa mappings
INSERT INTO aa_codes values('g' , 'gly');
INSERT INTO aa_codes values('a' , 'ala');
INSERT INTO aa_codes values('v' , 'val');
INSERT INTO aa_codes values('l' , 'leu');
INSERT INTO aa_codes values('i' , 'ile');
INSERT INTO aa_codes values('m' , 'met');
INSERT INTO aa_codes values('f' , 'phe');
INSERT INTO aa_codes values('w' , 'trp');
INSERT INTO aa_codes values('p' , 'pro');
INSERT INTO aa_codes values('s' , 'ser');
INSERT INTO aa_codes values('t' , 'thr');
INSERT INTO aa_codes values('c' , 'cys');
INSERT INTO aa_codes values('y' , 'tyr');
INSERT INTO aa_codes values('n' , 'asn');
INSERT INTO aa_codes values('q' , 'gln');
INSERT INTO aa_codes values('d' , 'asp');
INSERT INTO aa_codes values('e' , 'glu');
INSERT INTO aa_codes values('k' , 'lys');
INSERT INTO aa_codes values('r' , 'arg');
INSERT INTO aa_codes values('h' , 'his');
INSERT INTO aa_codes values('1' , 'pca');

-- delete existing data
TRUNCATE pgd_core_protein;
TRUNCATE pgd_core_residue;
TRUNCATE pgd_core_chain;

-- import proteins
INSERT INTO pgd_core_protein
    SELECT
        code        AS code,
        threshold   AS threshold,
        resolution  AS resolution,
        rfactor     AS rfactor
    FROM protein_info;

-- import chains
INSERT INTO pgd_core_chain (id, protein_id, code)
    SELECT DISTINCT 
        CONCAT(code, chainID) AS id,
        code                  AS protein_id,
        chainID               AS code
    FROM protein;

-- import residues
INSERT INTO pgd_core_residue (protein_id,chain_id,aa,chainID,oldID,chainIndex,a1,a2,a3,a4,a5,a6,a7,L1,L2,L3,L4,L5,ss,phi,psi,chi,ome,bm,bs,bg,h_bond_energy,zeta,terminal_flag,xpr)
    SELECT
        code            AS protein_id,
        CONCAT(code, chainID)   AS chain_id,
        aa_codes.key    AS aa,
        chainID         AS chainID,
        oldID           AS oldID,
        id              AS chainIndex,
        a1              AS a1,
        a2              AS a2,
        a3              AS a3,
        a4              AS a4,
        a5              AS a5,
        a6              AS a6,
        a7              AS a7,
        L1              AS L1,
        L2              AS L2,
        L3              AS L3,
        L4              AS L4,
        L5              AS L5,
        ss              AS ss,
        phi             AS phi,
        psi             AS psi,
        chi             AS chi,
        ome             AS ome,
        bm              AS bm,
        bs              AS bs,
        bg              AS bg,
        H_bond_energy   AS h_bond_energy,
        Zeta            AS zeta,
        terminal_flag   AS terminal_flag,
        xpr             AS xpr

    FROM protein, aa_codes
    where protein.AA = upper(aa_codes.value)

