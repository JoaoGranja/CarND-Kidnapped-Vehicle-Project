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
  num_particles = 10;  // TODO: Set the number of particles
  
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  std::cout << "GRANJA Init " << "x = " << x << " y =  " << y << " theta = " << theta << std::endl; 
  
  for (int i = 0; i< num_particles; i++)
  {
    Particle particle;
    
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.theta = remainder(particle.theta, (2 * M_PI));
    particle.weight = 1;  
    
    particles.push_back(particle);
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
  
  for (int i = 0; i< num_particles; i++)
  {   
    
    double theta     = particles[i].theta;
    double new_theta = theta + yaw_rate*delta_t; 
    particles[i].x = particles[i].x + ((velocity/yaw_rate) * (sin(new_theta) - sin(theta)) + dist_x(gen));
    particles[i].y = particles[i].y + ((velocity/yaw_rate) * (cos(theta) - cos(new_theta)) + dist_y(gen));
    particles[i].theta = new_theta + dist_theta(gen);
    particles[i].theta = remainder(particles[i].theta, (2 * M_PI)); 
    
    std::cout << "GRANJA predict " << "x = " << particles[i].x << " y =  " << particles[i].y << " theta = " << particles[i].theta << std::endl; 
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
 
  
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
       	observations[j].id = predicted[i].id;
      }
      else if (min_distance > distance )
      {
         min_distance = distance;
         observations[j].id = predicted[i].id;
      }
    }
  }

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
  
  for (int i = 0; i < num_particles; i++)
  { 
    
    // Transformation
    vector<LandmarkObs> trans_observations;
  	for (unsigned int j = 0; j < observations.size(); j++)
  	{
    	LandmarkObs obs;
    	obs.x = particles[i].x + (cos(particles[i].theta)*observations[j].x) - (sin(particles[i].theta)*observations[j].y);
    	obs.y = particles[i].y + (sin(particles[i].theta)*observations[j].y) + (cos(particles[i].theta)*observations[j].x);
      	trans_observations.push_back(obs);
      
      std::cout << "GRANJA trans_observations " << "x = " << obs.x << " y =  " << obs.y << std::endl;
    }
       
    //std::cout << "GRANJA TESTE before association" << std::endl;
    
    // Association
    vector<LandmarkObs> landmarks_obs;
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
  	{
      LandmarkObs obs;
      obs.x = map_landmarks.landmark_list[j].x_f;
      obs.y = map_landmarks.landmark_list[j].y_f;
      obs.id = map_landmarks.landmark_list[j].id_i - 1;
      landmarks_obs.push_back(obs);
      
      std::cout << "GRANJA landmarks_obs " << "x = " << obs.x << " y =  " << obs.y << std::endl;
    }

    dataAssociation(landmarks_obs, trans_observations);
    
    double prob = 1.0;
    for (unsigned int k = 0; k < trans_observations.size(); k++)
  	{
      // mult-variate Gaussian distribution
      //std::cout << "GRANJA TESTE before dist" << particles[i].x << particles[i].y << transformed_observations[k].x << transformed_observations[k].y << std::endl;
      double distance = dist(particles[i].x, particles[i].y, trans_observations[k].x, trans_observations[k].y); 
      //std::cout << "GRANJA TESTE after dist" << std::endl;
      
      if (distance <= sensor_range)
      {
        double x_obs = trans_observations[k].x;
        double y_obs = trans_observations[k].y;
        double mu_x  = landmarks_obs[trans_observations[k].id].x;
  		double mu_y  = landmarks_obs[trans_observations[k].id].y;
        
        //std::cout << "GRANJA TESTE before multiv_prob" << std::endl;
    	prob *= multiv_prob(std_landmark[0], std_landmark[1],  x_obs, y_obs, mu_x, mu_y);
        std::cout << "GRANJA multiv_prob " << prob << " " << x_obs << " " << mu_x << " " << y_obs << " " << mu_y<< std::endl;  
        //std::cout << "GRANJA TESTE after multiv_prob" << std::endl;
      }
      //std::cout << "GRANJA TESTE aqui" << std::endl;
    }
    //std::cout << "GRANJA TESTE weight " << prob << std::endl;
    particles[i].weight = prob;
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
  
  std::vector<Particle> resample_particles;

  auto max_it = std::max_element(particles.begin(), particles.end(), 
                                 [] (const Particle& p1, const Particle& p2) {return p1.weight < p2.weight;});
  double max_w = max_it->weight;
  int index = rand() % num_particles;
  double B = 0;
  
  for (int i = 0; i < num_particles; i++)
  {
      B = B + ((2*max_w) * rand());
      while(particles[index].weight < B)
      {
          B = B - particles[index].weight;
          index = remainder((index + 1),  num_particles);
      }
      resample_particles.push_back(particles[index]);
  }

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