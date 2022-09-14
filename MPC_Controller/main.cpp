/*!
 * @file main.cpp
 * @brief Main Function for the WBC Controller
 *
 * The main function parses command line arguments and starts the appropriate
 * driver.
 */

#include <main_helper.h>
#include "Imitation_Controller.hpp"

int main(int argc, char** argv) {
  main_helper(argc, argv, new Imitation_Controller());
  return 0;
}
