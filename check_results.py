import pickle as pkl 
import plotext as plt
import sys

job_name = sys.argv[1]

with open(f'./{job_name}/loc_min_history.pkl', 'rb') as h:
    loc_history = pkl.load(h)

with open(f'./{job_name}/all_history.pkl', 'rb') as h:
    all_history = pkl.load(h)

rmsd_all = [el[0] for el in all_history.values()]
rmsd_loc = [el[0] for el in loc_history.values()]
p0, rmsd0 = next(iter(all_history.items()))
print(f'default parameters: {p0}')
print(f'rmsd_0 = {rmsd0[0]}')
rmsd_combined = [rmsd_all[0]] + rmsd_loc
print(f'The number of all calculations: {len(rmsd_all)}')
'''
x_all = range(len(rmsd_all))
plt.plot(x_all, rmsd_all, marker='braille')
plt.xticks(x_all)
plt.plotsize(50, 15)
plt.title("all historical values")
plt.show()
plt.clear_figure()
'''
x_combined = range(len(rmsd_combined))
plt.plot(x_combined, rmsd_combined, marker='braille')
plt.xticks(x_combined)
plt.plotsize(50, 15)
plt.title("local minimas")
plt.show()

print(rmsd_combined)
last_key, last_value = list(loc_history.items())[-1]
print()
print('best solution:', last_key)
print('the best value:', last_value[0])
