import pyconf
import random
import tensorflow as tf

# -- configuration --

maxT = 1000
minT = 500
dT   = 10

num_hidden_units = 80

num_iterations = 100
num_equilibrate_steps_per_site = 100
num_sample_steps_per_site = 1

# -------------------

if __name__ == '__main__':

  N = len(pyconf.getSpins())
  num_y = len(range(minT, maxT, dT))

  x  = tf.placeholder(tf.float32, shape=[None, N])
  y_ = tf.placeholder(tf.float32, shape=[None, num_y])

  W_h = tf.Variable(tf.zeros([N, num_hidden_units]))
  W_o = tf.Variable(tf.zeros([num_hidden_units, num_y]))

  b_h = tf.Variable(tf.zeros([num_hidden_units]))
  b_o = tf.Variable(tf.zeros([num_y]))

  # Create the model
  def model(X, w_h, b_h, w_o, b_o):
      h = tf.sigmoid(tf.matmul(X, w_h) + b_h)
      pyx = tf.nn.softmax(tf.matmul(h, w_o) + b_o)
      return pyx

  y_hypo = model(x, W_h, b_h, W_o, b_o)

  # Cost Function basic term
  cross_entropy = -tf.reduce_sum(y_*tf.log(y_hypo))

  # Regularization terms (weight decay)
  L2_sqr = tf.nn.l2_loss(W_h) + tf.nn.l2_loss(W_o)
  lambda_2 = 0.01

  # the loss and accuracy
  loss = cross_entropy + lambda_2 * L2_sqr
  train_step = tf.train.GradientDescentOptimizer(0.1).minimize(loss)
  correct_prediction = tf.equal(tf.argmax(y_hypo,1), tf.argmax(y_,1))
  accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))

  sess = tf.InteractiveSession()
  tf.global_variables_initializer().run()


  for iteration_step  in range( 1, num_iterations+1 ):

    sampled_labels = []
    sampled_spins  = []

    for index, T in enumerate( range(minT, maxT, dT) ):

      arr_index = [ 0 for _ in range(minT, maxT, dT) ]
      arr_index[index] = 1

      pyconf.setEquilibriate(T, N*num_equilibrate_steps_per_site)

      # sample
      for _ in range(N*num_sample_steps_per_site):
        pyconf.metropolis_step(T)
        sampled_spins.append(pyconf.getSpins())
        sampled_labels.append(arr_index)

    c = list(zip(sampled_spins, sampled_labels))
    random.shuffle(c)
    sampled_spins, sampled_labels = zip(*c)

    print(
      str(iteration_step) + " step : " +
      "accuracy = " +
      str(sess.run(accuracy, feed_dict={x: sampled_spins , y_: sampled_labels}))
    )

    for i in range(0, N-1):
      print( len(sampled_spins)/N )
      sess.run(train_step, feed_dict={
        x:  sampled_spins[i*len(sampled_spins)/N : (i+1)*len(sampled_spins)/N],
        y_: sampled_labels[i*len(sampled_labels)/N : (i+1)*len(sampled_labels)/N]
      })


  f = open("weight.out", "w")
  for row, w in enumerate(sess.run(W_o)) :
    for col, i in enumerate(w) :
      f.write( str(col)+" "+str(row)+" "+str(i)+"\n" )
    f.write("\n")
  f.close()
