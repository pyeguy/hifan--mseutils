'''
sqlite3 implementation of storage class for mse specs

'''

import sqlite3
from . import *
import json


CREATE_TBL_SQL = """
    CREATE TABLE {tbl_name}(
    idx INTEGER PRIMARY KEY,
    sampid TEXT,
    rt REAL,
    ccs REAL,
    mz REAL,
    ppm REAL,
    z INTEGER,
    n INTEGER,
    i REAL,
    ms2_data TEXT,
    mgf_files TEXT,
    src_frag_ids TEXT
    )

    """
CREATE_INDEX_SQL = """
    CREATE INDEX feature_idx ON {tbl_name}({col_name});
    """
INSERT_SQL = """
    INSERT INTO {name}(
        idx,
        sampid,
        rt,
        ccs,
        mz,
        ppm,
        z,
        n,
        i,
        ms2_data,
        mgf_files,
        src_frag_ids
        )
    VALUES(
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?,
        ?
        )
    """


def create_mse_table(conn,tbl_name):
    try:
        cur = conn.cursor()
        cur.execute(CREATE_TBL_SQL.format(tbl_name=tbl_name))
        index_cols = "sampid,rt,ccs,mz"
        cur.execute(CREATE_INDEX_SQL.format(tbl_name=tbl_name,col_name=index_cols))
    except Exception as e:
        print(e)


def create_db(dbname):
    conn = sqlite3.connect("{dbname}.db".format(dbname=dbname))
    with conn:
        create_mse_table(conn,'mse_specs')
        create_mse_table(conn,'source_frags')
    return conn

def _make_val_tup(idx,mse):
    vals = (
        idx,
        mse.sampid,
        mse.rt.val,
        mse.ccs.val,
        mse.mz.mz,
        mse.mz.ppm,
        mse.mz.z,
        len(mse.mgf_files), #this shouldn't be necissary..
        mse.i,
        json.dumps([[mz.mz,i] for mz,i in mse.ms2_data.items()]),
        json.dumps(list(mse.mgf_files)), #coersce to list for json ser
        json.dumps([idx+i+1 for i in range(len(mse.src_frags))])
        )
    return idx,vals

def add_mse(cur,mse,idx,sampid):
    idx,vals = _make_val_tup(idx,mse)
    cur.execute(INSERT_SQL.format(name="mse_specs"),vals)
    for srcfrg in mse.src_frags:
        srcfrg.sampid = sampid
        idx += 1
        idx,vals = _make_val_tup(idx,srcfrg)
        cur.execute(INSERT_SQL.format(name="source_frags"),vals)
    return idx

def add_mses(conn,mses,idx=0,sampid='n/a'):
    with conn:
        cur = conn.cursor()
        # non-atomic commits for (huge) perf boost
        cur.execute("BEGIN")
        for mse in mses:
            mse.sampid = sampid
            idx = 1 + add_mse(cur,mse,idx,sampid=sampid)
        cur.execute("COMMIT")
    return idx


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def gen_mses(dbname,tbl_name='mse_specs',query=None):
    conn = sqlite3.connect(dbname)
    with conn:
        conn.row_factory = dict_factory
        cur = conn.cursor()
        if query:
            cur = cur.execute(query)
        else:
            cur = cur.execute("SELECT * FROM {tbl_name}".format(tbl_name=tbl_name))
        pbar = tqdm(total=cur.rowcount)
        frag_cur =conn.cursor()
        while True:
            row_dicts = cur.fetchmany(1000)
            if row_dicts:
                for rd in row_dicts:
                    mse = MseSpec.from_sqlite_dict(rd)
                    if mse.src_frag_ids:
                        frag_ids = '('+ ",".join(map(str,mse.src_frag_ids)) + ")"
                        frag_cur.execute(
                            """SELECT * FROM source_frags 
                            WHERE idx IN {frag_ids}"""
                            .format(frag_ids=frag_ids)
                            )
                        for frag_row in frag_cur.fetchall():
                            mse.src_frags.append(MseSpec.from_sqlite_dict(frag_row))
                    pbar.update()
                    yield mse
            else:
                break


