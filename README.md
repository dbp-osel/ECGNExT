# ECGNExT
General Information

This software code outputs an estimate of a noise signal collected as part of an electrocardiograph (ECG) recording. The objective is to obtain representative samples of noise and motion artifacts from ECG devices under test (DUT) that can then be used for testing the robustness of ECG analysis algorithms. The algorithm for noise estimation is described in [1]. The core operation of this software is based on removing the instances of the QRS complex from the ECG so that an estimate of the noise component of the recording remains. Required inputs are the ECG signal recorded with noise/motion artifacts and locations (indices) of the R-peaks.

Please refer to the "ECGnoiseExtraction_Instructions_for_use_V01.pdf" file for further details and citation information.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Update: 06/27/2024

The ECG data used to analyze the performance of this algorithm is available in the "ECG_files" directory of this      repository along with the document describing the data (ECG_data_description.pdf). The analysis results are published in [2].

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Regulatory Science Tool (RST) Reference
•	RST Reference Number: RST24IP01.01  
•	Date of Publication: 08/08/2023  
•	Recommended Citation: U.S. Food and Drug Administration. (2023). ECG Noise Extraction Tool (ECGNExT) (RST24IP01.01). https://cdrh-rst.fda.gov/ecg-noise-extraction-tool-ecgnext  

# Disclaimer
About the Catalog of Regulatory Science Tools  
The enclosed tool is part of the Catalog of Regulatory Science Tools, which provides a peer- reviewed resource for stakeholders to use where standards and qualified Medical Device Development Tools (MDDTs) do not yet exist. These tools do not replace FDA-recognized standards or MDDTs. This catalog collates a variety of regulatory science tools that the FDA's Center for Devices and Radiological Health's (CDRH) Office of Science and Engineering Labs (OSEL) developed. These tools use the most innovative science to support medical device development and patient access to safe and effective medical devices. If you are considering using a tool from this catalog in your marketing submissions, note that these tools have not been qualified as Medical Device Development Tools and the FDA has not evaluated the suitability of these tools within any specific context of use. You may request feedback or meetings for medical device submissions as part of the Q-Submission Program.  

For more information about the Catalog of Regulatory Science Tools, email OSEL_CDRH@fda.hhs.gov.


This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

[1] L. Galeotti and C. G. Scully, "A method to extract realistic artifacts from electrocardiogram recordings for robust algorithm testing," Journal of Electrocardiology, vol. 51, no. 6, Supplement, pp. S56-S60, 2018/11/01/ 2018, doi: 10.1016/j.jelectrocard.2018.08.023.

[2] A. Suliman, M. Farahmand, L. Galeotti, and C. G. Scully, "Clinical Evaluation of the ECG Noise Extraction Tool as a Component of ECG Analysis Algorithms Evaluation," IEEE Transactions on Biomedical Engineering, pp. 1-10, 2024, doi: 10.1109/TBME.2024.3386493.
