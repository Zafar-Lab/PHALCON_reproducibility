import pickle
with open('posMap.pkl', 'rb') as f:
    loaded_dict = pickle.load(f)
print(loaded_dict)