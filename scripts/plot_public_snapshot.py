"""Regenerate figures that require only shipped aggregate outputs."""
import plot_phase3_results as p
def main():
    p.setup_style()
    evaluations=p.pd.read_csv(p.ROOT/'results/performance/phase1_v1/study_specific_evaluations.csv')
    evaluations=evaluations[evaluations.execution_status=='EXECUTED'].copy()
    tasks=p.pd.read_csv(p.ROOT/'metadata/phase1/phase1_task_summary.csv')
    tasks=tasks[tasks.status=='EXECUTABLE'].copy()
    inventory=p.pd.read_csv(p.ROOT/'literature/signature_inventory.csv')
    inventory=inventory[inventory.main_grid_eligible.astype(str).str.upper()=='TRUE'].copy()
    cells=p.pd.read_csv(p.RESULTS/'phase3_cell_robustness.csv')
    p.plot_figure_1(evaluations,tasks,inventory);p.plot_figure_3(evaluations);p.plot_figure_4(cells);p.plot_figure_5();p.plot_figure_6(cells);p.plot_supplementary(evaluations,cells)
    print('Rebuilt Figures 1, 3–6 and S1–S4. Figure 2 distribution panel needs locally regenerated subject scores.')
if __name__=='__main__':main()
