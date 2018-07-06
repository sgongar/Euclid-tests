#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
from numpy import array
import pandas as pd
from sklearn import preprocessing, svm
from sklearn.preprocessing import StandardScaler
from sklearn import model_selection, metrics
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split, KFold, cross_val_score, cross_val_predict

import matplotlib.pyplot as plt
from matplotlib import style
import datetime

style.use('ggplot')
raw_data = 'catalogue.csv'
df = pd.read_csv(raw_data, index_col=0)

# Create a DataFrame for numerical features
# data1 = pd.DataFrame(df, columns=['PM', 'A_IMAGE', 'B_IMAGE', 'CLASS_STAR'])

# Create a DataFrame for categorical features
# cols_to_transform = pd.DataFrame(df, columns=['Forest', 'Closed_Shrublands', 'Open_Shrublands',
#                                              'Woody_Savannas', 'Savannas', 'Grasslands', 'Croplands'])
# dummies = pd.get_dummies(cols_to_transform)
# Join data1 and dummies using Numpy and yield as array
X = array(pd.DataFrame(df, columns = ['B_IMAGE', 'A_IMAGE']))

# Specify the dependent variable as array
y = array(df['PM'])


lm = LinearRegression(n_jobs=-1)

"""
To check the accuracy / confidence level of the prediction,
we have 25 % test datasets, while 75 % is used for training.
"""

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25)
print (X_train.shape, y_train.shape)
print (X_test.shape, y_test.shape)
# First we fit a model
model = lm.fit(X_train, y_train)
# print the coefficents
print("The linear cofficients", model.coef_)
# Try to predict the y ( NPP_Predict) for the test data-features(independent variables(X_test)
predictions = lm.predict(X_test)
# Accuracy of the prediction
confidence = lm.score(X_test, y_test)
print("This is predicted NPP2001 Values", predictions)
print("This is the prediction accuracy", confidence)

plt.legend(loc=4)
plt.title("Actual NPP2001 vs. NPP2001_Predict", size = 10)
plt.scatter(y_test, predictions, color='c', marker ='.')
plt.xlabel("Actual NPP2001", size = 10)
plt.ylabel("NPP2001_Predict", size = 10)
plt.show()

plt.legend(loc=4)
plt.title("Homogeneity of Variance")
plt.scatter(y_test, y_test - predictions)
plt.xlabel("Actual NPP2001")
plt.ylabel("Residual")
plt.show()
# Perform 10 fold Cross Validation (KFold)
scores = cross_val_score(model, X, y, cv=10)
print ("Cross Validated Scores", scores)
kf = KFold(n_splits=10, random_state=None, shuffle=True)
for train_index, test_index in kf.split(X):
    print ("TRAIN", train_index, "TEST", test_index)
X_train, X_test = X[train_index], X[test_index]
y_train, y_test = y[train_index], y[test_index]
# Make Cross Validated predictions
predictions2 = cross_val_predict(model, X, y, cv=10)
# Check the R2- the proportion of variance in the dependent variable explained by the predictors
accuracy = metrics.r2_score(y, predictions2)
print("This is R2", accuracy)
plt.scatter(y, predictions2, color='c', marker ='.')
plt.legend(loc=4)
plt.xlabel("Actual NPP2001", size = 10)
plt.ylabel("NPP2001_Predict", size = 10)
plt.title("Actual and Predicted NPP2001 Values using 10 Fold Cross Validation", size = 10)
plt.show()
