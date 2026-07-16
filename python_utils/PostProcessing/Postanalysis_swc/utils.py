
import os

def Export_folders_from_results(PATH):
    """Detect all subfolders inside PATH (simulation results)."""
    simulation_dirs = sorted([
        os.path.join(PATH, d)
        for d in os.listdir(PATH)
        if os.path.isdir(os.path.join(PATH, d))
    ])
    for s in simulation_dirs:
        print(s)
    return simulation_dirs
