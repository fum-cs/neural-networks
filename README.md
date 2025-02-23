![](lectures/img/572_banner.png)

## Computer Science Dept., Ferdowsi University of Mashhad

# Neural Networks

An introduction to Neural Networks with Python and Pytorch which covers optmization, neural network basics, convolutional neural networks, and advanced topics such as autoencoders and generative adversarial networks.

- [Course Jupyter Book](https://fum-cs.github.io/neural-networks/README.html)

2024 Instructor: Mahmood Amintoosi

---

I should mention that the original material was from [Tomas Beuzen's course github](https://github.com/UBC-MDS/DSCI_572_sup-learn-2). I have modified his contents to suit my own needs and preferences. I would like to thank him for his great work and generosity.


## Optional Reference/Learning Materials

### Deep learning resources
- [Dive into Deep Learning](http://d2l.ai/chapter_introduction/index.html), a book based on STAT 157 at UC Berkeley.
- [Deep learning YouTube series](https://www.youtube.com/watch?v=aircAruvnKk) by 3Blue1Brown.
- [Neural Networks and Deep Learning](http://neuralnetworksanddeeplearning.com/) (free online book).
- [Deep Learning](http://www.deeplearningbook.org/). Ian Goodfellow, Yoshua Bengio and Aaron Courville.
- [Deep Learning with Python](https://machinelearningmastery.com/deep-learning-with-python). Jason Brownlee.
- [Stanford UFLDL tutorial](http://deeplearning.stanford.edu/wiki/index.php/UFLDL_Tutorial) (or [here](http://deeplearning.stanford.edu/tutorial/))
- [Geoff Hinton Coursera lectures](https://www.youtube.com/playlist?list=PLoRl3Ht4JOcdU872GhiYWf6jwrk_SNhz9)
- [CS231n: Convolutional Neural Networks for Visual Recognition (Stanford)](http://cs231n.github.io/)
- [Grokking Deep Learning](https://www.manning.com/books/grokking-deep-learning)
- [Practical Deep Learning For Coders, Part 1](http://course.fast.ai/) and some more resources on their blog [here](http://www.fast.ai/2016/12/19/favorite-posts/)
- [A Guide to Deep Learning](http://yerevann.com/a-guide-to-deep-learning/)
- [Awesome Deep Learning](https://github.com/ChristosChristofidis/awesome-deep-learning), which is a list of other resources
- [Full Stack Deep Learning](https://fullstackdeeplearning.com/)

### ML-related textbooks

- Hands-On Machine Learning with Scikit-Learn, Keras, and TensorFlow: Concepts, Tools, and Techniques to Build Intelligent Systems by Aurélien Géron. Code/notebooks available [here](https://github.com/ageron/handson-ml2). (Endorsed by an MDS student!)
- James, Gareth; Witten, Daniela; Hastie, Trevor; and Tibshirani, Robert. An Introduction to Statistical Learning: with Applications in R. 2014. Plus [Python code](https://github.com/JWarmenhoven/ISLR-python) and [more Python code](https://github.com/mscaudill/IntroStatLearn).
- Russell, Stuart, and Peter Norvig. Artificial intelligence: a modern approach. 1995.
- David Poole and Alan Mackwordth. Artificial Intelligence: foundations of computational agents. 2nd edition (2017). [Free e-book](http://artint.info/).
- Kevin Murphy. Machine Learning: A Probabilistic Perspective. 2012.
- Christopher Bishop. Pattern Recognition and Machine Learning. 2007.
- Pang-Ning Tan, Michael Steinbach, Vipin Kumar. Introduction to Data Mining. 2005.
- Mining of Massive Datasets. Jure Leskovec, Anand Rajaraman, Jeffrey David Ullman. 2nd ed, 2014.

### Math for ML

- [Mathematics for Machine Learning](https://mml-book.github.io/)
- [The Matrix Calculus You Need For Deep Learning](http://parrt.cs.usfca.edu/doc/matrix-calculus/index.html)
- [Introduction to Optimizers](https://blog.algorithmia.com/introduction-to-optimizers/)

### Other ML resources

- [A Course in Machine Learning](http://ciml.info/)
- [Nando de Freitas lecture videos](https://www.youtube.com/watch?v=PlhFWT7vAEw) and [online course](https://www.cs.ox.ac.uk/people/nando.defreitas/machinelearning/)

### Interesting ML Competition Write-Ups

- [Diabetic retinopathy Kaggle competition write-up](http://jeffreydf.github.io/diabetic-retinopathy-detection/)
- [Galaxy Zoo Kaggle competition write-up](https://benanne.github.io/2014/04/05/galaxy-zoo.html)
- [National Data Science Bowl competition write-up](https://benanne.github.io/2015/03/17/plankton.html)


## Build

In notebooks folder:
- jupyter-book build ./
- copy ../require.js ./_build
- ghp-import -n -p -f ./_build/html
- jupyter-book build --builder pdflatex ./
