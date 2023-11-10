# Weekly Meeting Notes

**October 20, 2023**

- Colon cancer exists in middle of pipe
    - Loaded with bacteria
- Some colon tissue & colon cancer tissue covered in bacteria
- Big question: Relationship between bacteria touching colon cancer & colon cancer itself?
- Looking carefully at cancerous tissue
    - To see if there is any bacteria
- Epithelial cells: Cells those lines tissue
    - For any organ that secrets things
        - Breast, prostate
    - Most cancers come from epithelial cells
    - Separates us from microbes inside of us
        - Special, thin barrier
            - Ability to absorb nutrients, excrete waste, keep out pathogens
- In colon cancer & in some normal tissue, some bacteria invade or get very close to epithelial cells
    - Biofilm
        - Sugary protein that is goo that helps bacteria adhere to surface
        - Help protect bacteria from antibiotics
        - Resistant to own immune system
        - More common in right side of colon
        - Ask question: Which came first? Bacteria nestle up next to epithelial cells to form cancer or cancer formed first and bacteria got attracted
- C. difficile was a part of these biofilms
    - Bacteria that has some properties
        - Forms spores released in stool, carried onto next host (ex: via surfaces)
        - Spores in healthy GI tract do not do anything
            - In sick patients, will vegetate and produce toxins
                - TcdA, TcdB
                    - Mediate disease, kill epithelial cells
    - Patients who have C. difficile infection
        - Cause is these toxins
        - Become susceptible while taking antibiotics for something else
            - Can go into gut and kill healthy bacteria
    - Was not previously identified was related to causing colon cancer
        - Can accelerate colon cancer in genetically susceptible mice
- Question: How does C. difficile do that?
    - If toxin is what mediate disease in infectious diarrhea, is toxin what in a different patient lead to cancer?
        - Wanted to test that where can remove several variables
            -How big is person, diet, bathroom frequency
- Test it out in dishes
    - Matrigel
    - Expose them just to toxin, taking out many variables from C. difficile
        - Is it possible for C. difficile toxin itself to effect epithelial cells to induce cancer?
- 3728T
    - Identification of sample patient, T for tumor
- Consortium from 30 bacteria derived from original slurry from patient 3728 or mice who got slurry 3728
- At step where have raw sequencing data
    - Each of cDNA got sequenced in both directions
    - During reverse transcription single strand mRNA becomes into double strand DNA
    - Parameters built into it
- 2 separate files for samples
    - Samples numbered 1-30
    - Each have R1 and R2 file
        - Annotated by project ID, initials of person, sample number, R1 & R2 for separate files (fastQ files, compressed with gzip)
- Storage & Memory
    - Memory
    - Align these 1 by 1, computationally intense, require way more memory
    - Use package: kallisto 
- Next Steps: 
    - Look through documentation of kallisto
    - Need to set up private GitHub
        - Find out how much data can be uploaded

**October 27, 2023**
- Showed Dr. Markham GitHub Repo and tutorials created for kallisto
- Next Steps: 
    - Look at 'sleuth' package in R
    - Move raw data onto external hard drive for access
    - Experiment importing & reading in data to Kallisto
        - Follow manual to do alignments & mapping

**November 3, 2023**
- No meeting as data is currently being transferred to external hard drive
