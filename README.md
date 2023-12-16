How human-AI feedback loops alter human perceptual, emotional and social judgements
Moshe Glickman1, 2 & Tali Sharot1, 2, 3
1 Affective Brain Lab, Department of Experimental Psychology, University College London, London, UK
2 Max Planck UCL Centre for Computational Psychiatry and Ageing Research, University College London, London, UK
3 Department of Brain and Cognitive Sciences, Massachusetts Institute of Technology, Cambridge, MA, USA
Correspondence authors: mosheglickman345@gmail.com, t.sharot@ucl.ac.uk

Artificial intelligence (AI) technologies are rapidly advancing, enhancing human capabilities across various domains spanning from finance to medicine. Despite their numerous advantages, AI systems can exhibit biases in judgments ranging from perception to emotion. Here, in a series of experiments (N=1,201), we reveal a feedback loop where human-AI interactions alter processes underlying human perceptual, emotional and social judgements, subsequently amplifying biases in humans. This amplification is significantly greater than observed in interactions between humans, due both to the tendency of AI systems to amplify biases and to how humans perceive AI systems. Participants are often unaware of the extent of the AI’s influence, rendering them more susceptible to it.  These findings reveal a mechanism wherein AI systems amplify human biases, which are further internalized by humans during human-AI interactions, triggering a snowball effect where small errors in judgment escalate into much larger ones.
Experiment 1
Background
Experiment 1 aimed to test how interactions with AI influence human decision-making. To this end, we first collected data in an emotion aggregation task in which participants are presented with an array of 12 faces, and classify the mean emotion expressed by them as 'more sad' or 'more happy'. Humans show a slight bias to respond ‘more sad’. An AI (convolutional neural network) trained on this slightly biased dataset, further amplified the bias. Next, we show that humans interacting with the biased AI, became more biased over time. This bias amplification does not occur in the network human-human network. 
Data Files:
Level 1
•	Exp1-Level1.csv
•	Exp1-Level1-SingleFace.csv
Level 2
•	Exp1-Level2-H-H-H.csv
Level 3
•	Exp1-Level3-Baseline.csv (No interaction)
•	Exp1-Level3-H-AI-H.csv (Human-AI interaction)
•	Exp1-Level3-H-AIH-H.csv (Human-AI perceived-as-human interaction)
•	Exp1-Level3-H-HAI-H.csv (Human-human perceived-as-AI interaction)
•	Exp1-Level3-H-H-H.csv (Human-human interaction)
Main analyses
•	Analyzing the tendency of participants to classify arrays as ‘more sad’ compared to a 0.5 baseline (Level 1)
•	Analyzing the tendency of participants to classify arrays as ‘more sad’ compared to a 0.5 baseline (Level 2)
•	Assessing the percentage of change in decisions across the conditions in Level 3.
•	Evaluating how the interaction with AI and humans affected participants bias across time.
Experiment 2
Background 
Experiment 2 sought to mimic a situation in which humans are not a-priori biased, but rather AI bias emerges for other reasons (e.g., under-represented classes). To this end, we employed a variant of the Random Dot Kinematogram task. Participants interacted with three three simple algorithms which provided either a systematically biased response (biased AI), an accurate response (accurate AI) or a noisy response (noisy AI). The results indicated that interaction with a biased AI resulted in a significant human bias relative to baseline and that interaction with an accurate AI resulted in a significant increase in accuracy relative to baseline. A reinforcement learning model provided a good fit for the data, and outperformed a baseline model which did not assume learning from the AI.
Data Files:
Exp2.csv – Interaction with the biased, accurate and noisy algorithms.
Exp2-Accurate.csv - Interaction with the accurate algorithm.
Exp2-Biased.csv - Interaction with the biased algorithm.
Main analyses
•	Calculating bias and accuracy for different conditions (Accurate AI, Biased AI, Noisy AI).
•	Analyzing participants' perceived influence of AI in different conditions.
•	A computational reinforcement learning model.
Experiment 3
Background
Experiment 3 tested whether exposure to an AI system with subtle gender biases can bias human social judgments. Participants evaluated pairs of applicants for a job requiring a specific ability. They then interacted with a biased AI model favoring men in its judgments of mixed-gender pairs. The results indicate that, as compared to women, participants evaluated men as more competent after interacting with a gender-biased AI. 
Data Files:
•	allDataRelative.csv – Relative judgments
•	allDataSeparate.csv – Separate judgments
Main analyses:
•	Analyzes of AI-induced biases in relative judgments
•	Analyzes of AI-induced biases in separate judgments
Experiment 4
Background:
Experiment 4 demonstrated a feedback loop using the outputs of a widely used real world text-to-image generative AI – Stable Diffusion. We first show that Stable Diffusion (trained on human-generated images) generates biased outputs. In particular, when asked to generate images of financial managers Stable Diffusion generates images of white man 85% of the time, which is a rate that exceeds their actual representation in the population (for example, U.S. – 44.3% men, 78.5% white overall; U.S. Bureau of Labor Statistics,2022; U.K. – 55.6% men; Office for National Statistics, 2021). We then show that interacting with Stable Diffusion’s biased outputs further biases participants’ judgements about the likelihood of white-men financial managers, beyond their initial judgements.
Data Files:
•	Exp4.csv
Main analyses:
•	Analyzing the tendency of the treatment and control groups to choose white men after exposure to the AI images or to the control images, in comparison to their preferences before exposure.
•	Examining changes in participants' judgements (before vs. after exposure to the images) using a mixed-model multinomial logistic regression (conducted using SPSS).
