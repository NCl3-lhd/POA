"""

Author:zhai1xiao
Create Date:2021-12-06
Update Date:2025-03-31
Modified by: wym6912

"""
import re
import argparse
import multiprocessing



def read_fasta(file_path):
    """
    读取fasta文件
    :param file_path: 路径
    :return: identifiers [] 序列id的列表
             sequences [] 序列的列表
    """
    identifiers = []
    sequences = []
    with open(file_path, 'r') as fin:
        curr = ""
        flag = False
        for line in fin:
            if line[0] == '>':
                if flag:
                    sequences.append(curr.upper())
                curr = ""
                flag = True
                identifiers.append(line[1:].strip())
            else:
                curr += line.strip()
        sequences.append(curr.upper())
    return identifiers, sequences


def preprocess(sequences):
    """
    处理fasta序列中的非法字符（非字母用N代替）

    :param sequences 待处理的序列列表
    :return: processed [] 处理后的序列列表
    """
    processed = []
    for i in sequences:
        tmp = re.sub("[^ACTGNUactgnu-]", "N", str(i).replace('U', 'T'))
        processed.append(tmp)
    return processed

def score_of(curr_clm: dict, matchS, mismatchS, gap1S, gap2S):
    """
    序列单列打分

    Args:
        sequences  待打分的序列
        matchS     match(碱基相同且非N)的得分
        mismatchS  mismatch(碱基非空非N，但不相同)的得分
        gap1S      gap比对碱基，N比对N，N比对碱基的得分
        gap2S      gap比对gap，gap比对N的得分

    Return:
        该列的sp score
    """
    match = 0  # nongap==nongap
    mismatch = 0  # nongap!=nongap
    gap1 = 0  # gap-nongap
    gap2 = 0
    A, C, G, T, N, dash = curr_clm['A'], curr_clm['C'], curr_clm['G'], curr_clm['T'], curr_clm['N'], curr_clm['-']
    match = (A * (A - 1) + C * (C - 1) + G * (G - 1) + T * (T - 1)) // 2
    mismatch = ((A + C) * (G + T)) + A * C + G * T
    gap1 = (A + C + G + T + N) * dash
    gap2 = (dash * (dash - 1) // 2) + ((N * (N - 1)) // 2) + (A + C + G + T) * N
    score = (match * matchS) + (mismatch * mismatchS) + (gap1 * gap1S) + (gap2 * gap2S)
    return score

def every_place_calc(lock, score_dict, place_start, place_end, seqs):
    '''
    每个位置计算

    Args:
        lock        进程锁
        score_dict  结果写入的字典
        place_start 起始位置
        place_end   终止位置
        seqs        序列总数
    Return:
        当前区域的 sp 值

    '''
    _score = 0.0
    for place in range(place_start, place_end):
        curr_clm: dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, '-': 0}
        for i in range(seqs):
            curr_clm[sequences[i][place]] += 1
        _score += score_of(curr_clm, matchScore, mismatchScore, gap1Score, gap2Score)
    with lock:
        score_dict['score'] += _score
        print("Calculate range [{}, {}) done.".format(place_start, place_end))

def evaluate(sequences, matchS, mismatchS, gap1S, gap2S):
    """
    序列打分

    Args:
        sequences  待打分的序列
        matchS     match(碱基相同)的得分
        mismatchS  mismatch(碱基非空，但不相同)的得分
        gap1S      gap比对碱基的得分
        gap2S      gap比对gap的得分

    Return:
        sp score
    """
    print("calculating", flush=True)
    score = 0.0
    '''
    for j in range(len(sequences[0])):
        curr_clm: dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, '-': 0}
        for i in range(len(sequences)):
            curr_clm[sequences[i][j]] += 1
        score += score_of(curr_clm, matchS, mismatchS, gap1S, gap2S)
        if (j + 1) % (len(sequences[0]) // 3) == 0: print('.', end="", flush=True)
    '''
    max_len = len(sequences[0])
    if(workers != 1):
        per = max_len // (workers - 1)
        rg = [per * i for i in range(workers)]
        rg.append(max_len)
        seqs = len(sequences)
        lock = multiprocessing.Lock()
        with multiprocessing.Manager() as score:
            score_dic = score.dict({'score': 0.0})
            results = [multiprocessing.Process(target = every_place_calc, args = (lock, score_dic, rg[i], rg[i + 1], seqs, )) for i in range(len(rg) - 1)]
            for _ in results:
                _.start()
            for _ in results:
                _.join()
            score = score_dic['score']
    else:
        per = max_len
        for j in range(len(sequences[0])):
            curr_clm: dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, '-': 0}
            for i in range(len(sequences)):
                curr_clm[sequences[i][j]] += 1
            score += score_of(curr_clm, matchS, mismatchS, gap1S, gap2S)
            if (j + 1) % (len(sequences[0]) // 3) == 0: print('.', end="", flush=True)
    print("done")
    return score


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=False, default='output.aligned.fasta', help='Fasta file path to be scored.')
    parser.add_argument('--match', type=float, default=1, help='match score, default=1.')
    parser.add_argument('--mismatch', type=float, default=-1, help='mismatch score, default=-1.')
    parser.add_argument('--gap1', type=float, default=-2, help='gap-base score, default=-2.')
    parser.add_argument('--gap2', type=float, default=0, help='gap-gap score, default=0.')
    parser.add_argument('--threads', type=int, default=1, help='program threads, default=1.')
    args = parser.parse_args()

    global matchScore
    global mismatchScore
    global gap1Score
    global gap2Score
    global sequences
    global workers

    filepath = args.input
    matchScore = args.match
    mismatchScore = args.mismatch
    gap1Score = args.gap1
    gap2Score = args.gap2
    workers = args.threads


    _, sequences = read_fasta(filepath)
    print(str(len(sequences)) + " sequences found")
    if len(sequences) <= 1: exit()

    sequences = preprocess(sequences)


    sp = evaluate(sequences, matchScore, mismatchScore, gap1Score, gap2Score)  # pure sp
    avgSp = sp / (len(sequences) * (len(sequences) - 1) // 2)  # avg sp
    scaledSP = avgSp / len(sequences[0])  # scaled sp
    print("SP score: " + str(sp))
    print("Avg SP score: " + str(avgSp))
    print("Scaled SP score: " + str(scaledSP))
    print("Finish!")