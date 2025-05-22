import os
import sys
import importlib.util
import matplotlib.pyplot as plt

def import_module_from_path(module_name, file_path):
    """Import a module from a file path."""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

def run_example(example_path, example_name):
    """Run an example script and handle its execution."""
    print(f"\n{'='*80}")
    print(f"Running {example_name}...")
    print(f"{'='*80}")
    
    # Import the module
    module_name = os.path.basename(example_path).replace('.py', '')
    try:
        module = import_module_from_path(module_name, example_path)
        
        # Run the main function of the module
        if hasattr(module, 'main'):
            module.main()
        else:
            print(f"Error: {example_name} does not have a main() function.")
    except Exception as e:
        print(f"Error running {example_name}: {str(e)}")
    
    # Close any open plots to avoid interference between examples
    plt.close('all')

def main():
    """Run all examples in the project."""
    print("ModelFit Examples Runner")
    print("=======================\n")
    
    # Get the project root directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    examples_dir = os.path.join(script_dir, 'examples')
    
    # Define the examples to run
    examples = [
        {
            'path': os.path.join(examples_dir, 'biomass_growth_model.py'),
            'name': 'Monod Growth Model'
        },
        {
            'path': os.path.join(examples_dir, 'biomass_growth_death_model.py'),
            'name': 'Growth-Death Model'
        },
        {
            'path': os.path.join(examples_dir, 'dual_substrate_growth_model.py'),
            'name': 'Dual Substrate Growth Model'
        }
    ]
    
    # Run each example
    for example in examples:
        run_example(example['path'], example['name'])
    
    print("\nAll examples completed!")
    print("Results saved in the 'figures' directory.")

if __name__ == "__main__":
    main()
