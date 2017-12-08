import datetime as dt

import pandas as pd

def parse_ldsc_rg_log(fn):
    conv_month = {'': 0, 'Apr': 4, 'Aug': 8, 'Dec': 12, 'Feb': 2, 
                  'Jan': 1, 'Jul': 7, 'Jun': 6, 'Mar': 3, 
                  'May': 5, 'Nov': 11, 'Oct': 10, 'Sep': 9}
    with open(fn) as f:
        fcontents = f.read()
    lines = fcontents.split(69 * '*' + '\n')[-1].strip().split('\n')
    month, day, time, year = [x.split() for x in lines if x[0:10] == 'Beginning '][0][4:]
    hour, minute, second = time.split(':')
    begin = dt.datetime(int(year), int(conv_month[month]), int(day), int(hour), int(minute), int(second))

    month, day, time, year = [x.split() for x in lines if x[0:17] == 'Analysis finished'][0][4:]
    hour, minute, second = time.split(':')
    end = dt.datetime(int(year), int(conv_month[month]), int(day), int(hour), int(minute), int(second))

    num_snps = int([x for x in lines if 'valid' in x][0].split()[0])

    # Pheno 1
    lines = fcontents.split(69 * '*' + '\n')[-1].split(29 * '-' + '\n')[0].strip().split('\n')
    p1_h2, p1_h2_se = [x for x in lines if x[0:5] == 'Total'][0].split()[-2:]
    p1_h2 = float(p1_h2)
    p1_h2_se = float(p1_h2_se[1:-1])
    p1_lambda_gc = float([x for x in lines if x[0:6] == 'Lambda'][0].strip().split()[-1])
    p1_mean_chi2 = float([x for x in lines if x[0:4] == 'Mean'][0].strip().split()[-1])
    p1_intercept, p1_intercept_se = [x for x in lines if x[0:9] == 'Intercept'][0].strip().split()[-2:]
    p1_intercept = float(p1_intercept)
    p1_intercept_se = float(p1_intercept_se[1:-1])
    
    # Pheno 2
    lines = fcontents.split(69 * '*' + '\n')[-1].split(29 * '-' + '\n')[0].strip().split('\n')
    p2_h2, p2_h2_se = [x for x in lines if x[0:5] == 'Total'][0].split()[-2:]
    p2_h2 = float(p2_h2)
    p2_h2_se = float(p2_h2_se[1:-1])
    p2_lambda_gc = float([x for x in lines if x[0:6] == 'Lambda'][0].strip().split()[-1])
    p2_mean_chi2 = float([x for x in lines if x[0:4] == 'Mean'][0].strip().split()[-1])
    p2_intercept, p2_intercept_se = [x for x in lines if x[0:9] == 'Intercept'][0].strip().split()[-2:]
    p2_intercept = float(p2_intercept)
    p2_intercept_se = float(p2_intercept_se[1:-1])
    
    vals = [begin, end, num_snps]
    ind = ['start_time', 'end_time', 'num_snps']
    vals += [p1_h2, p1_h2_se, p1_lambda_gc, p1_mean_chi2, p1_intercept, 
             p1_intercept_se]
    ind += ['h2_p1', 'h2_se_p1', 'lambda_gc_p1', 'mean_chi2_p1', 'intercept_p1',
           'intercept_se_p1']
    vals += [p2_h2, p2_h2_se, p2_lambda_gc, p2_mean_chi2, p2_intercept, 
             p2_intercept_se]
    ind += ['h2_p2', 'h2_se_p2', 'lambda_gc_p2', 'mean_chi2_p2', 'intercept_p2',
           'intercept_se_p2']
    lines = fcontents.split(69 * '*' + '\n')[-1].strip().split('\n')
    vals += lines[-4].split()[0:2]
    ind += lines[-5].split()[0:2]
    vals += [float(x) for x in lines[-4].split()[2:]]
    ind += lines[-5].split()[2:]
    
    out = pd.Series(vals, index=ind)
    return(out)
