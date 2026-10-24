# BEAD 当前仓库架构

> 基于 2026-10-22 工作区中的实际代码生成。实线表示主要运行时调用或数据流，虚线表示配置、资源或状态读取。

## 总体架构图

```mermaid
flowchart TB
    U["用户 / Python 调用"] --> APP["bead.app.BEAD<br/>命令行交互与单次调用入口"]

    subgraph ORCH["Agent 编排层 · bead/agents"]
        AGENT["BEADAgent"]
        PROMPT["bead.prompt<br/>核心提示词 + 动态任务上下文"]
        GRAPH["LangGraph StateGraph<br/>MessagesState"]
        MODELNODE["agent 节点<br/>ChatOpenAI adapter"]
        TOOLNODE["tools 节点<br/>ToolNode + tools_condition"]
        MEMORY["MemorySaver<br/>按 thread_id 保存进程内对话状态"]

        AGENT --> GRAPH
        PROMPT -. SystemMessage .-> MODELNODE
        GRAPH --> MODELNODE
        MODELNODE -->|"tool_calls"| TOOLNODE
        TOOLNODE -->|"ToolMessage / observation"| MODELNODE
        MEMORY <-.-> GRAPH
    end

    APP --> AGENT
    MODELNODE <-->|"OpenAI API · HTTPS"| OPENAI["OpenAI GPT-4.1 API<br/>API key 来自 OPENAI_API_KEY"]
    MODELNODE -->|"最终回复"| APP

    TOOLNODE --> REGISTRY["bead.agents.make_tools<br/>16 个 LangChain 工具包装器"]

    subgraph TOOLS["领域工具实现层 · bead/tools"]
        CHEM["tools.py<br/>子结构检索 / 反向过滤<br/>3D 相似性 / 组合工作流"]
        ADMET["tools.py · ADMET_prediction<br/>本地 ADMET 推理与校准"]
        MOLPREP["molecule_preparation.py<br/>小分子标准化与 3D 准备"]
        TARGET["target_preparation.py<br/>蛋白准备与对接盒生成"]
        GEN["generation.py + vina_rl.py<br/>REINVENT / LibINVENT 强化学习"]
        SDOCK["docking.py<br/>Schrödinger LigPrep + Glide"]
        ADOCK["autodock_docking.py<br/>AutoDock Vina 对接"]
        STATUS["5 类 job_status 工具<br/>读取 status.json / job.log"]
    end

    REGISTRY --> CHEM
    REGISTRY --> ADMET
    REGISTRY --> MOLPREP
    REGISTRY --> TARGET
    REGISTRY --> GEN
    REGISTRY --> SDOCK
    REGISTRY --> ADOCK
    REGISTRY --> STATUS

    subgraph LOCAL["仓库内数据与模型资源"]
        MOLLIST["bead/source/mollist.csv<br/>分子-靶点活性库"]
        SCAFFOLD["bead/source/Scaffold_1.mol"]
        ADMETRES["bead/resources/admet_ai<br/>ADMET-AI checkpoints<br/>半衰期 ensemble<br/>CL priority calibrator<br/>DrugBank reference"]
        GENRES["bead/resources/generation<br/>REINVENT / LibINVENT priors<br/>TOML 模板"]
    end

    CHEM --> MOLLIST
    CHEM --> SCAFFOLD
    CHEM --> RDKIT["RDKit + bead.utils.sim.Sim3D"]
    ADMET --> ADMETRES
    ADMET --> ADMETAI["ADMET-AI 2.0.1 + RDKit/joblib"]
    MOLPREP --> RDKIT
    GEN --> GENRES

    subgraph BACKENDS["计算后端"]
        SCH["Schrödinger 2024-1<br/>Protein Prep / LigPrep / Glide"]
        OPEN["开源准备栈<br/>PDBFixer / OpenMM<br/>Molscrub / Meeko"]
        VINA["AutoDock Vina<br/>ProLIF / MDAnalysis"]
        RL["Maize workflow<br/>REINVENT4 / LibINVENT"]
        GPU["PyTorch / CUDA GPU"]
    end

    TARGET -->|"优先"| SCH
    TARGET -->|"自动回退"| OPEN
    SDOCK --> SCH
    ADOCK --> OPEN
    ADOCK --> VINA
    GEN --> RL
    RL -->|"后端 A"| SCH
    RL -->|"后端 B"| VINA
    RL --> GPU

    TARGET -->|"prepared PDB/PDBQT<br/>Glide grid / Vina box"| SDOCK
    TARGET -->|"prepared PDB/PDBQT<br/>Vina box"| ADOCK
    TARGET -->|"grid 或 receptor + box"| GEN
    MOLPREP -->|"SDF/PDB/PDBQT"| ADOCK

    subgraph OUTPUTS["输出与可观测性"]
        SYNC["同步结果<br/>CSV / PNG"]
        JOB["后台任务目录<br/>job.json · status.json · job.log"]
        ART["计算产物<br/>CSV · SDF · PDB/PDBQT<br/>交互报告/图片 · checkpoints"]
    end

    CHEM --> SYNC
    ADMET --> SYNC
    MOLPREP --> JOB
    TARGET --> JOB
    GEN --> JOB
    SDOCK --> JOB
    ADOCK --> JOB
    JOB --> ART
    STATUS -.轮询.-> JOB
```

## Agent 的 ReAct 运行流程

```mermaid
sequenceDiagram
    actor User as 用户
    participant App as bead.app.BEAD
    participant Agent as BEADAgent / LangGraph
    participant LLM as OpenAI GPT-4.1 API
    participant Tool as ToolNode / BEAD tools
    participant Job as 后台 worker

    User->>App: 自然语言任务 + thread_id
    App->>Agent: chat()/run()
    Agent->>Agent: 组装系统提示词与 MessagesState
    Agent->>LLM: SystemMessage + 对话历史
    LLM-->>Agent: 回复或 tool_calls
    alt 需要调用同步工具
        Agent->>Tool: 参数化工具调用
        Tool-->>Agent: CSV 路径 / 计算结果
    else 需要耗时任务
        Agent->>Tool: 启动 generation/docking/preparation
        Tool->>Job: 独立 subprocess + job.json
        Tool-->>Agent: job_dir / status_file / pid
        loop 状态查询
            Agent->>Tool: *_job_status_tool(job_path)
            Tool-->>Agent: 状态、日志尾部、产物列表
        end
    end
    Agent->>LLM: ToolMessage / observation
    LLM-->>Agent: 下一步工具调用或最终答案
    Agent-->>App: 最终自然语言回复
    App-->>User: 结果与文件路径
```

## 活跃工具映射

| Agent 工具 | 主要实现 | 核心输入 | 主要输出/后端 |
|---|---|---|---|
| `Subsearch` / `Rev_subsearch` | `bead/tools/tools.py` | scaffold、CSV | RDKit 子结构过滤 CSV/图片 |
| `Similarity_prediction` | `tools.py` + `utils/sim.py` | reference ID | RDKit ETKDG/MMFF 3D 相似性 CSV |
| `filter_then_similarity_tool` | `tools.py` | ID、scaffold、相似性范围 | 靶点过滤 → 反向过滤 → 3D 相似性 |
| `ADMET_prediction` | `tools.py` | SMILES CSV | 本地 ADMET-AI + 专用半衰期模型 + CL 校准器 |
| `Mol_preparation` | `molecule_preparation.py` | SMILES/InChI/CSV | RDKit/Meeko；SDF、PDB、PDBQT |
| `Target_preparation` | `target_preparation.py` | PDB、可选 box | Schrödinger 或 PDBFixer/OpenMM/Meeko |
| `Denovo_molgen` | `generation.py` | grid 或 receptor/box | REINVENT + Glide/Vina RL；默认 100 epochs |
| `Constrained_molgen` | `generation.py` | scaffold + grid 或 receptor/box | LibINVENT + Glide/Vina RL；默认 1500 epochs |
| `Docking_schrodinger` | `docking.py` | SMILES CSV + Glide grid | LigPrep/Glide、排序与相互作用报告 |
| `Docking_autodock` | `autodock_docking.py` | SMILES CSV + receptor + box | Molscrub/Meeko/Vina/ProLIF |
| 5 个 `*_job_status_tool` | 对应后台模块 | job 目录或 `status.json` | 状态、日志尾部、产物路径 |

## 仓库目录职责

```text
BEAD/
├── bead/
│   ├── app.py                    # 对外入口
│   ├── prompt.py                 # Agent 系统提示词
│   ├── agents/
│   │   ├── agent.py              # LangGraph ReAct 编排、OpenAI GPT-4.1 API 连接
│   │   └── make_tools.py         # LangChain 工具注册与参数接口
│   ├── tools/                    # 化学计算、生成、准备、对接及状态查询
│   ├── utils/sim.py              # 3D 相似性底层实现
│   ├── resources/                # 本地模型、训练后校准器、生成先验
│   └── source/                   # 默认分子库与 scaffold
├── scripts/                      # 数据准备、半衰期/CL 训练与分析脚本
├── vendor/reinvent4/             # 仓库内 REINVENT4 源码副本
├── tests/                        # 测试及 1018 实验运行数据
├── benchmarks/                   # Agent 性能评测
├── setup.py / environment.yml    # Python 3.11 包与环境定义
└── *_TOOL.md                     # 各重型工具的使用文档
```

## 当前实现边界

- 默认大模型是 OpenAI API 的 `gpt-4.1`，凭据从 `OPENAI_API_KEY` 环境变量读取。
- `MemorySaver` 是进程内存记忆；进程退出后不会持久化对话状态。
- 生成、对接、分子准备和靶点准备采用独立后台进程；Agent 通过状态工具读取日志和产物。
- `bead/tools/tools.py` 是当前 Agent 使用的检索、相似性和 ADMET 主实现；`tools h.py`、`tools h2.py`、`tools2.py` 等未注册到当前 `bead_toolset`，属于历史/实验文件。
- Schrödinger 是可选闭源后端；准备和生成支持自动回退到开源栈，独立的 `Docking_schrodinger` 本身仍要求 Schrödinger。
