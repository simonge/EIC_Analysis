
import tensorflow.keras.backend
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras.optimizers import Adam

#tensorflow.keras.backend.set_floatx('float64')
 
# Define model
model = Sequential()
model.add(Dense(16,   activation='tanh', input_shape=(4,)))
model.add(Dense(4,    activation='tanh'))
model.add(Dense(1,    activation='linear'))
 
# Set loss and optimizer
model.compile(loss='mean_squared_error', optimizer=Adam(), weighted_metrics=[])
 
# Store model to file
model.save('modelRegression.h5')
model.summary()
