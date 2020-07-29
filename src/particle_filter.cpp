/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  
  for (int i = 0; i< num_particles; i++)
  {
    Particle particle;
    
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;  
    
    particles.push_back(particle);
    
    weights.push_back(particle.weight);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0.0, std_pos[0]);
  std::normal_distribution<double> dist_y(0.0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0.0, std_pos[2]);
  
  //std::cout << "prediction" << std::endl;
  
  for (int i = 0; i< num_particles; i++)
  {   
    // in case yaw rate is too small simplify the predictions equations
    if (fabs(yaw_rate) < 0.001) 
    {
      particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
      particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
      particles[i].theta += dist_theta(gen);
    }
    else
    {
      double theta     = particles[i].theta;
      double new_theta = theta + yaw_rate*delta_t;
      particles[i].x += ((velocity/yaw_rate) * (sin(new_theta) - sin(theta))) + dist_x(gen);
      particles[i].y += ((velocity/yaw_rate) * (cos(theta) - cos(new_theta))) + dist_y(gen);
      particles[i].theta = new_theta + dist_theta(gen);
    }
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations,
                                     Particle &particle) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */ 
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;
  
    //std::cout << "dataAssociation" << std::endl;
  
  for (unsigned int j = 0; j < observations.size(); j++)
  {
    bool first_time = true;
    long distance, min_distance;
    for (unsigned int i = 0; i < predicted.size(); i++)
  	{
      distance = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
      if (first_time == true)
      {
       	first_time = false;
       	min_distance = distance;
       	observations[j].id = i;
      }
      else if (min_distance > distance )
      {
         min_distance = distance;
         observations[j].id = i;
      }
    }
    associations.push_back(predicted[observations[j].id].id);
    sense_x.push_back(observations[j].x);
    sense_y.push_back(observations[j].y);
  }
  
  // Set assocations for visualization purpose only
  SetAssociations(particle, associations, sense_x, sense_y);

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  //std::cout << "updateWeights" << std::endl;
  
  double total_weight = 0;
  for (int i = 0; i < num_particles; i++)
  { 
    // Transformation of observation from car coordinates to map coordinates 
    vector<LandmarkObs> trans_observations;
  	for (unsigned int j = 0; j < observations.size(); j++)
  	{
    	LandmarkObs obs;
      	obs.id = observations[j].id;
    	obs.x = particles[i].x + (cos(particles[i].theta)*observations[j].x) - (sin(particles[i].theta)*observations[j].y);
    	obs.y = particles[i].y + (sin(particles[i].theta)*observations[j].x) + (cos(particles[i].theta)*observations[j].y);
      	trans_observations.push_back(obs);
    }
    
    // Getting map landmarks positions and id within the sensor range
  	vector<LandmarkObs> landmarks_obs;
  	for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
  	{
    	LandmarkObs obs;   
    	obs.x = map_landmarks.landmark_list[j].x_f;
    	obs.y = map_landmarks.landmark_list[j].y_f;
    	obs.id = map_landmarks.landmark_list[j].id_i;
      
    	double distance = dist(particles[i].x, particles[i].y, obs.x, obs.y);
      	//Compute within the sensor range
      	if (distance <= sensor_range)
      	{	
    		landmarks_obs.push_back(obs);
        }
  	}    
       
    // landmark association with observation
    dataAssociation(landmarks_obs, trans_observations, particles[i]);
    
    // mult-variate Gaussian distribution
    double prob = 1.0;
    
    for (unsigned int k = 0; k < trans_observations.size(); k++)
  	{
      double mu_x  = landmarks_obs[trans_observations[k].id].x;
  	  double mu_y  = landmarks_obs[trans_observations[k].id].y;

      double x_obs = trans_observations[k].x;
      double y_obs = trans_observations[k].y;
      double weight = multiv_prob(std_landmark[0], std_landmark[1],  x_obs, y_obs, mu_x, mu_y); 
        
      if(weight > 0)
      	prob *= weight;
    }
    
    particles[i].weight = prob;
    total_weight += prob;
  }
  
  
  //Weights normaliation
  for (int i = 0; i < num_particles; i++)
  { 
    particles[i].weight = particles[i].weight/ total_weight;
    weights[i] = particles[i].weight;
  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  /* initialize random seed: */
  srand (time(NULL));
  
  //std::cout << "resample" << std::endl;
  
  std::vector<Particle> resample_particles;

  auto max_it = std::max_element(particles.begin(), particles.end(), 
                                 [] (const Particle& p1, const Particle& p2) {return p1.weight < p2.weight;});
  double max_w = max_it->weight;
  
  int index = rand() % num_particles;
  double B = 0;
  
  for (int i = 0; i < num_particles; i++)
  {
      B = B + ((2*max_w) * ((double) rand() / (RAND_MAX)));
      while(B > particles[index].weight)
      {
          B = B - particles[index].weight;
          index = fmod((index + 1), num_particles);
      }
      resample_particles.push_back(particles[index]);
  }
  particles = resample_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}