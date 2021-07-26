# coding: utf-8
import uproot
import numpy as np
import os

signal = np.random.multivariate_normal([2.0, 5.0], [[1.0, 0.7], [0.7, 1.0]], size=1000000)
background = np.random.multivariate_normal([5.0, 2.0],[[1.0, -0.7],[-0.7, 1.0]], size=1000000)
sfile = uproot.recreate("assets/signal.root")
sfile["bat"] = uproot.newtree({"energy":float, "position":float}, title="Batman")
sfile["bat"].extend({"energy":signal.T[0], "position":signal.T[1]})
sfile.close()
bfile = uproot.recreate("assets/background.root")
bfile["bat"] = uproot.newtree({"energy":float, "position":float}, title="Batman")
bfile["bat"].extend({"energy":background.T[0], "position":background.T[1]})
bfile.close()
