# dynamic_report/ directory

This directory holds the subpipeline that takes care of producing the dynamic portable report document. Having an independent subpipeline that runs after the main "CompareM2" workflow, was the best way I could make something that was robust against intermittently failing jobs. 
